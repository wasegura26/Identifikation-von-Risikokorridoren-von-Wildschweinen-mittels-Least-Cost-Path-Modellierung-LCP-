#==============================================#
# Identifikation von Risikokorridoren von Wildschweinen mittels Least-Cost-Path LCP
# im Kontext der Prävention der Afrikanischen Schweinepest ASP – Fallstudie Landkreis Miltenberg
# Autor: Walter Segura – Nov 2025
#==============================================#

#------------------------------------------------
# 1. Initialisierung
# Ziel: Arbeitsumgebung bereinigen und Fachpakete laden
# Methode: Speicher und Grafiken leeren sowie S2 in sf deaktivieren
#------------------------------------------------
rm(list = ls()); if (!is.null(dev.list())) graphics.off(); cat("\014")
library(terra)
library(sf)
library(leastcostpath)
try(sf::sf_use_s2(FALSE), silent = TRUE)

# Projektwurzel sowie portables Temp-Verzeichnis für terra
projekt_root <- "D:/Project/Wildschwein_LCP_Miltenberg"
temp_dir     <- file.path(projekt_root, "Temp_Terra")
if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
terraOptions(tempdir = temp_dir, memfrac = 0.6, progress = TRUE)

#------------------------------------------------
# 2. Parameterdefinition
# Ziel: Zentrale Bündelung aller Pfade und Modellparameter
#------------------------------------------------
pfad_landbedeckung <- file.path(projekt_root, "Input", "Miltenberg_Raster_Produktlayer_CLMS_10m.tif")
pfad_hangneigung   <- file.path(projekt_root, "Input", "Miltenberg_Slope_10m.tif")
pfad_einstaende    <- file.path(projekt_root, "Input", "Zentroide_Sommer_Winter_Einstand.shp")

ordner_ausgabe <- file.path(projekt_root, "Output", "Ergebnisse_LCP")
if (!dir.exists(ordner_ausgabe)) dir.create(ordner_ausgabe, recursive = TRUE)

ziel_crs            <- "EPSG:25832"   # metrisches Projektionssystem
suchradius_m        <- 7000           # ökologisch motivierter Suchradius
anzahl_kandidaten   <- 200            # Begrenzung der Zielkandidaten zwecks Effizienz
hang_faktor         <- 0.08           # Sensitivität der Steigung im Exponentialterm
nachbarschaft       <- 8              # Nachbarschaft für die Konnektivität
faktor_abs          <- 1.75           # absolute Längenbegrenzung
faktor_rel          <- 2.00           # relative Längenbegrenzung
klasse_wasser       <- 10             # Zielklasse Wasser
klasse_nahrung      <- 12             # Zielklasse Acker
wasser_beta         <- 1.0            # Intensität der Wasserfolgenstrafe
na_klassen          <- c()            # zusätzliche Barrieren als NA
puffer_siedlung_m   <- 10             # Sperrpuffer um Siedlung in Meter

# Ausgabewege als Shapefile für unmittelbares Feedback
fmt_w <- function(ii) file.path(ordner_ausgabe, sprintf("Einstand_%02d_Wasserweg.shp",  ii))
fmt_n <- function(ii) file.path(ordner_ausgabe, sprintf("Einstand_%02d_Nahrungsweg.shp", ii))
pfad_alle_shp <- file.path(ordner_ausgabe, "Alle_Wege.shp")

#------------------------------------------------
# 3. Datenimport und Harmonisierung
# Ziel: Einheitliches Koordinatenreferenzsystem sowie identische Auflösung und Ausdehnung
# Methode: Projektion anpassen und Hangraster bilinear resamplen
#------------------------------------------------
r_land       <- rast(pfad_landbedeckung)  # Landbedeckung kategorisch
r_hang       <- rast(pfad_hangneigung)    # Hangneigung kontinuierlich in Grad
punkte_einst <- vect(pfad_einstaende)     # Einstände als Punkte

if (crs(r_land)       != ziel_crs) r_land       <- project(r_land,       ziel_crs)
if (crs(r_hang)       != ziel_crs) r_hang       <- project(r_hang,       ziel_crs)
if (crs(punkte_einst) != ziel_crs) punkte_einst <- project(punkte_einst, ziel_crs)

if (!isTRUE(all.equal(res(r_land), res(r_hang)))) r_hang <- resample(r_hang, r_land, method = "bilinear")
if (!isTRUE(all.equal(ext(r_land), ext(r_hang)))) { r_hang <- crop(r_hang, r_land); r_hang <- resample(r_hang, r_land, method = "bilinear") }

# Teilmenge der Einstände zur gezielten Prüfung
nur_auswahl   <- TRUE
index_bereich <- 31:31
index_einst   <- if (nur_auswahl) intersect(index_bereich, seq_len(nrow(punkte_einst))) else seq_len(nrow(punkte_einst))
if (!length(index_einst)) stop("Keine Einstände im gewählten Bereich gefunden.")

#------------------------------------------------
# 4. Kostenmodell
# Ziel: Zellenkosten definieren und harte Barrieren festlegen
# Formel: Kosten gleich Basis mal Exponentialterm der Steigung
#------------------------------------------------
reklass_matrix <- matrix(c(
  0,  0, NA,     # NoData gilt als Barriere
  1,  1, NA,     # Siedlung gilt als Barriere
  2,  2, 1,      # Nadelwald
  3,  3, 1,      # Laubwald
  5,  5, 1,      # Gehölz
  6,  6, 100,    # Dauergrünland
  7,  7, 20,     # Krautige Vegetation
  9,  9, 1000,   # Offenland
  10, 10, 5000,  # Wasser Basiswert für die Folgenstrafe
  12, 12, 1      # Acker
), ncol = 3, byrow = TRUE)

kosten_basis <- classify(r_land, rcl = reklass_matrix, include.lowest = TRUE)

# Siedlungspuffer als harte Barriere
if (puffer_siedlung_m > 0) {
  siedlung_bin <- ifel(r_land == 1, 1, NA)
  d2sied       <- distance(siedlung_bin)
  mask_s       <- !is.na(d2sied) & d2sied <= puffer_siedlung_m
  kosten_basis[mask_s] <- NA
}

# Optionale Zusatzbarrieren
if (length(na_klassen) > 0) kosten_basis[r_land %in% na_klassen] <- NA

# Gesamtkosten der Landbewegung mit Propagation von NA
kosten_raster        <- kosten_basis * exp(hang_faktor * r_hang)

# Strikte Maskierung auf begehbare Zellen gemäß Basis
kosten_raster_masked <- mask(kosten_raster, kosten_basis)

#------------------------------------------------
# 5. Konnektivität und Komponentenlogik
# Ziel: Konnektivitätsobjekt erzeugen und begehbare Komponenten kennzeichnen
# Methode: create_cs anwenden sowie zusammenhängende Patches identifizieren
#------------------------------------------------
konnektiv <- create_cs(kosten_raster_masked, neighbours = nachbarschaft)

walkable <- !is.na(kosten_raster_masked)  # wahr auf begehbaren Zellen
comp_r   <- terra::patches(walkable, directions = if (nachbarschaft %in% c(8,16,32)) 8 else 4)

get_comp_id <- function(pt, comp_r){
  v <- terra::extract(comp_r, pt, cells = FALSE)[,2]
  if (length(v) == 0 || is.na(v[1])) return(NA_integer_) else return(as.integer(v[1]))
}
subset_same_component <- function(pts, comp_r, comp_id){
  if (is.na(comp_id) || is.null(pts) || nrow(pts) == 0) return(NULL)
  vv <- terra::extract(comp_r, pts, cells = FALSE)[,2]
  keep <- which(!is.na(vv) & vv == comp_id)
  if (!length(keep)) return(NULL)
  pts[keep]
}

#------------------------------------------------
# 6. Zielableitung Wasser und Nahrung
# Ziel: Uferzellen als Wasserziele und Ackerzellen als Nahrungsziele bestimmen
# Methode: Ufer über lokale Nachbarschaft sowie Acker über Klassenwert ableiten
#------------------------------------------------
wmask <- r_land == klasse_wasser
nb    <- focal(wmask, w = matrix(1,3,3), fun = sum, na.rm = TRUE, pad = TRUE)
ufer  <- wmask & nb < 9 & nb > 0

if (global(ufer, "sum", na.rm = TRUE)[1,1] == 0) stop("Keine Uferzellen gefunden.")
ziel_w_xy     <- xyFromCell(ufer, which(values(ufer) == 1))
ziele_wasser  <- if (length(ziel_w_xy)) terra::vect(ziel_w_xy, type = "points", crs = crs(r_land)) else NULL

zellen_acker  <- which(values(r_land) == klasse_nahrung)
ziel_n_xy     <- xyFromCell(r_land, zellen_acker)
ziele_nahrung <- if (length(zellen_acker)) terra::vect(ziel_n_xy, type = "points", crs = crs(r_land)) else NULL
if (is.null(ziele_nahrung) || nrow(ziele_nahrung) == 0) stop("Keine Ackerzellen vorhanden.")

#------------------------------------------------
# 7. Hilfsfunktionen für Liniengeometrien
# Ziel: Gültigkeit sichern und minimale Attributierung vornehmen
#------------------------------------------------
linie_gueltig <- function(vlin){
  if (is.null(vlin) || nrow(vlin) == 0) return(FALSE)
  s  <- sf::st_as_sf(vlin)
  if (any(sf::st_is_empty(s))) return(FALSE)
  gt <- as.character(sf::st_geometry_type(s, by_geometry = TRUE))
  if (!any(gt %in% c("LINESTRING","MULTILINESTRING"))) return(FALSE)
  ln <- suppressWarnings(as.numeric(sf::st_length(s)))
  if (all(is.na(ln) | ln == 0)) return(FALSE)
  TRUE
}
annotiere_linien <- function(vlin, id_einstand, ziel_txt, ress_txt){
  g <- sf::st_as_sf(vlin)
  g$einstand_id <- as.integer(id_einstand)
  g$ziel        <- as.character(ziel_txt)
  g$ressource   <- as.character(ress_txt)
  terra::vect(g)
}

#------------------------------------------------
# 8. Pfadkostenfunktion mit Wasserfolgenlogik
# Ziel: Gesamtkosten eines Pfades bestimmen
# Methode: Summe der Landkosten sowie exponentielle Zuschläge bei aufeinanderfolgenden Wasserzellen
# Randbedingung: NA in Basis oder Kosten gilt als Barriere
#------------------------------------------------
pfad_kosten_expo_streak_fast <- function(vlin, r_land, r_kosten, r_basis, klasse_wasser, wasser_beta = 1.0){
  if (is.null(vlin) || nrow(vlin) == 0) return(Inf)
  if (!identical(terra::crs(vlin), terra::crs(r_land))) vlin <- terra::project(vlin, terra::crs(r_land))
  
  extc <- try(terra::extract(r_land, vlin, cells = TRUE), silent = TRUE)
  if (inherits(extc, "try-error") || is.null(extc) || nrow(extc) == 0) return(Inf)
  
  cells <- unique(extc$cell)
  lv    <- values(r_land)[cells]
  kv    <- values(r_kosten)[cells]
  bv    <- values(r_basis)[cells]
  
  if (any(is.na(bv))) return(Inf)
  if (any(is.na(kv[lv != klasse_wasser]))) return(Inf)
  
  xy       <- xyFromCell(r_land, cells)
  cs       <- res(r_land)[1]
  step_len <- rep(1, length(cells))
  if (length(cells) > 1){
    dxy <- sqrt(rowSums((xy[-1,,drop=FALSE] - xy[-nrow(xy),,drop=FALSE])^2))
    step_len[-1] <- dxy / cs
  }
  
  total <- 0; racha <- 0
  for (i in seq_along(cells)){
    if (lv[i] == klasse_wasser){
      racha  <- racha + 1
      total  <- total + bv[i] * exp(wasser_beta * racha) * step_len[i]
    } else {
      racha  <- 0
      total  <- total + kv[i] * step_len[i]
    }
  }
  total
}

#------------------------------------------------
# 9. Auswahlregel für den besten Pfad
# Ziel: Für jeden Einstand den Pfad mit minimalen Gesamtkosten zu Wasser und Acker bestimmen
# Methode: Filter auf begehbare Komponente sowie Radius und Top K dann LCP dann Guardrails dann Kostenbewertung
#------------------------------------------------
bester_pfad <- function(startpunkt, zielpunkte, konn, radius_m, top_k,
                        faktor_abs, faktor_rel,
                        r_land, r_kosten, r_basis, klasse_wasser, wasser_beta = 1.0){
  comp_id_start <- get_comp_id(startpunkt, comp_r)
  if (is.na(comp_id_start)) {
    message("  · Ursprung außerhalb der begehbaren Zone, Eintrag wird übersprungen.")
    return(NULL)
  }
  zielpunkte <- subset_same_component(zielpunkte, comp_r, comp_id_start)
  if (is.null(zielpunkte) || nrow(zielpunkte) == 0) {
    message("  · Keine Ziele im gleichen Komponentenverbund, Eintrag wird übersprungen.")
    return(NULL)
  }
  d_alle <- terra::distance(startpunkt, zielpunkte)
  idx_in <- which(is.finite(d_alle) & d_alle <= radius_m)
  if (!length(idx_in)) return(NULL)
  idx_sort <- idx_in[order(d_alle[idx_in])]
  idx_wahl <- idx_sort[1:min(top_k, length(idx_sort))]
  
  best_line <- NULL; best_cost <- Inf
  for (j in seq_along(idx_wahl)){
    ziel_j <- zielpunkte[idx_wahl[j]]
    if (is.na(get_comp_id(ziel_j, comp_r))) next
    
    pfad_j <- tryCatch(
      create_lcp(konn, origin = startpunkt, destination = ziel_j,
                 cost_distance = FALSE, check_locations = TRUE),
      error = function(e) NULL
    )
    if (is.null(pfad_j) || nrow(pfad_j) == 0) next
    
    laenge_j <- tryCatch(sum(as.numeric(sf::st_length(sf::st_as_sf(pfad_j))), na.rm = TRUE),
                         error = function(e) Inf)
    d_j    <- as.numeric(d_alle[idx_wahl[j]])
    grenze <- max(faktor_abs * radius_m, faktor_rel * d_j)
    if (!(is.finite(laenge_j) && laenge_j <= grenze)) next
    
    c_j <- pfad_kosten_expo_streak_fast(
      vlin = pfad_j, r_land = r_land, r_kosten = kosten_raster_masked, r_basis = kosten_basis,
      klasse_wasser = klasse_wasser, wasser_beta = wasser_beta
    )
    if (is.finite(c_j) && c_j < best_cost){
      best_cost <- c_j; best_line <- pfad_j
    }
  }
  if (!is.null(best_line) && linie_gueltig(best_line)) return(best_line)
  NULL
}

#------------------------------------------------
# 10. Hauptberechnung und inkrementelle Ausgabe
# Ziel: Für jeden Einstand getrennte Pfade Wasser und Acker berechnen und sofort schreiben
#------------------------------------------------
alle_pfade <- list()
for (ii in index_einst){
  start_i <- punkte_einst[ii]
  message(sprintf("[Einstand %02d] Suche nach besten Pfaden läuft …", ii))
  
  pfad_w <- bester_pfad(
    startpunkt = start_i, zielpunkte = ziele_wasser, konn = konnektiv,
    radius_m = suchradius_m, top_k = anzahl_kandidaten,
    faktor_abs = faktor_abs, faktor_rel = faktor_rel,
    r_land = r_land, r_kosten = kosten_raster_masked, r_basis = kosten_basis,
    klasse_wasser = klasse_wasser, wasser_beta = wasser_beta
  )
  if (linie_gueltig(pfad_w)){
    pfad_w <- annotiere_linien(pfad_w, ii, "wasser", "wasser")
    terra::writeVector(pfad_w, fmt_w(ii), filetype = "ESRI Shapefile", overwrite = TRUE)
    alle_pfade[[length(alle_pfade) + 1]] <- pfad_w
    message(sprintf("[Einstand %02d] Wasserweg geschrieben.", ii))
  } else {
    message(sprintf("[Einstand %02d] Kein gültiger Wasserweg.", ii))
  }
  
  pfad_n <- bester_pfad(
    startpunkt = start_i, zielpunkte = ziele_nahrung, konn = konnektiv,
    radius_m = suchradius_m, top_k = anzahl_kandidaten,
    faktor_abs = faktor_abs, faktor_rel = faktor_rel,
    r_land = r_land, r_kosten = kosten_raster_masked, r_basis = kosten_basis,
    klasse_wasser = klasse_wasser, wasser_beta = wasser_beta
  )
  if (linie_gueltig(pfad_n)){
    pfad_n <- annotiere_linien(pfad_n, ii, "nahrung", "acker")
    terra::writeVector(pfad_n, fmt_n(ii), filetype = "ESRI Shapefile", overwrite = TRUE)
    alle_pfade[[length(alle_pfade) + 1]] <- pfad_n
    message(sprintf("[Einstand %02d] Nahrungsweg geschrieben.", ii))
  } else {
    message(sprintf("[Einstand %02d] Kein gültiger Nahrungsweg.", ii))
  }
}

#------------------------------------------------
# 11. Sammelausgabe
# Ziel: Alle gültigen Pfade in einer Datei zusammenführen
# Hinweis: Einzeldateien bleiben zur Kontrolle bestehen
#------------------------------------------------
if (length(alle_pfade) > 0){
  alle_pfade_sv <- tryCatch(do.call(rbind, alle_pfade), error = function(e) NULL)
  if (!is.null(alle_pfade_sv) && inherits(alle_pfade_sv, "SpatVector")){
    terra::writeVector(alle_pfade_sv, pfad_alle_shp, filetype = "ESRI Shapefile", overwrite = TRUE)
    message("Sammelausgabe geschrieben: ", pfad_alle_shp)
  }
}

message("Fertig.")
