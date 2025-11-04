#!/usr/bin/env python3
import json
import os
import sys
import hashlib
from pathlib import Path
from typing import Dict, Iterable, List, Tuple, Optional
from urllib.parse import urljoin

import requests
import geopandas as gpd
from shapely.geometry import LineString, MultiLineString, Polygon, MultiPolygon, GeometryCollection
from shapely.ops import unary_union
import shapely
import pyproj
import mercantile
import mapbox_vector_tile

# -------------------------------
# User config
# -------------------------------

SKI_ROUTES_PATH = os.environ.get("SKI_ROUTES_PATH", "ski_routes.geojson")
Z = int(os.environ.get("MVT_ZOOM", "13"))

DANGER_LAYER_ALLOWLIST = os.environ.get("DANGER_LAYERS", "")
if DANGER_LAYER_ALLOWLIST.strip():
    DANGER_LAYER_ALLOWLIST = [s.strip() for s in DANGER_LAYER_ALLOWLIST.split(",")]
else:
    DANGER_LAYER_ALLOWLIST = []

PROPERTY_FILTER_JSON = os.environ.get("PROPERTY_FILTER_JSON", "").strip()
PROPERTY_FILTER: Optional[Dict[str, List[str]]] = json.loads(PROPERTY_FILTER_JSON) if PROPERTY_FILTER_JSON else None

THRESH_RED_MIN_CONTINUOUS = float(os.environ.get("RED_MIN_CONTINUOUS_M", "40"))
THRESH_GREEN_MAX_CONTINUOUS = float(os.environ.get("GREEN_MAX_CONTINUOUS_M", "12"))

OUTPUT_GEOJSON = os.environ.get("OUTPUT_GEOJSON", "ski_routes_classified_multi.geojson")

CACHE_DIR = Path(os.environ.get("MVT_CACHE_DIR", ".mvt_cache"))
CACHE_DIR.mkdir(parents=True, exist_ok=True)

READ_FROM_MBTILES = {
    "low": os.environ.get("READ_FROM_MBTILES_LOW", "").strip(),
    "moderate": os.environ.get("READ_FROM_MBTILES_MODERATE", "").strip(),
    "scary_moderate": os.environ.get("READ_FROM_MBTILES_SCARY_MODERATE", "").strip(),
    "considerable": os.environ.get("READ_FROM_MBTILES_CONSIDERABLE", "").strip(),
    "high": os.environ.get("READ_FROM_MBTILES_HIGH", "").strip(),
}

DEFAULT_DANGER_TILESETS = {
    "low": "https://api.maptiler.com/tiles/019756ec-8856-70ea-b772-9eb267efcf51/tiles.json?key=8vJgZ3VXNGNEoiXkOJxM",
    "moderate": "https://api.maptiler.com/tiles/019756fc-11b3-771e-b5a7-fa979f04813d/tiles.json?key=8vJgZ3VXNGNEoiXkOJxM",
    "scary_moderate": "https://api.maptiler.com/tiles/0197570d-e344-7419-8b9c-8d6d948bed3e/tiles.json?key=8vJgZ3VXNGNEoiXkOJxM",
    "considerable": "https://api.maptiler.com/tiles/01975725-3746-7a32-9bc0-ca3375eed5dc/tiles.json?key=8vJgZ3VXNGNEoiXkOJxM",
    "high": "https://api.maptiler.com/tiles/01975a28-a8bd-7f50-a3cc-2d53e43783ea/tiles.json?key=8vJgZ3VXNGNEoiXkOJxM",
}
DANGER_TILESET_JSON = os.environ.get("DANGER_TILESET_JSON", "").strip()
if DANGER_TILESET_JSON:
    try:
        DEFAULT_DANGER_TILESETS = json.loads(DANGER_TILESET_JSON)
    except Exception as e:
        print(f"[WARN] Could not parse DANGER_TILESET_JSON; using defaults. Error: {e}")

# -------------------------------
# Helpers
# -------------------------------

WEBM = pyproj.CRS.from_epsg(3857)
WGS84 = pyproj.CRS.from_epsg(4326)
PROJECT_TO_WEBM = pyproj.Transformer.from_crs(WGS84, WEBM, always_xy=True).transform

def tile_bounds_webm(z: int, x: int, y: int):
    b = mercantile.xy_bounds(mercantile.Tile(x=x, y=y, z=z))
    return (b.left, b.bottom, b.right, b.top)

def mvt_geom_to_webm(z: int, x: int, y: int, geom: Dict, extent: int = 4096):
    minx, miny, maxx, maxy = tile_bounds_webm(z, x, y)
    sx = (maxx - minx) / extent
    sy = (maxy - miny) / extent

    def _convert_coords(coords):
        return [(minx + sx * cx, maxy - sy * cy) for (cx, cy) in coords]

    gtype = geom["type"]
    if gtype == "Point":
        return shapely.Point(_convert_coords([geom["coordinates"]])[0])
    elif gtype == "MultiPoint":
        return shapely.MultiPoint(_convert_coords(geom["coordinates"]))
    elif gtype == "LineString":
        return shapely.LineString(_convert_coords(geom["coordinates"]))
    elif gtype == "MultiLineString":
        return shapely.MultiLineString([_convert_coords(line) for line in geom["coordinates"]])
    elif gtype == "Polygon":
        rings = [ _convert_coords(ring) for ring in geom["coordinates"] ]
        if not rings:
            return None
        return shapely.Polygon(rings[0], rings[1:] if len(rings) > 1 else None)
    elif gtype == "MultiPolygon":
        polygons = []
        for poly in geom["coordinates"]:
            rings = [ _convert_coords(ring) for ring in poly ]
            if rings:
                polygons.append(shapely.Polygon(rings[0], rings[1:] if len(rings) > 1 else None))
        return shapely.MultiPolygon(polygons) if polygons else None
    else:
        return None

def resolve_pbf_template(url: str) -> str:
    """
    Accept either a direct {z}/{x}/{y}.pbf template or a tiles.json URL.
    If tiles.json, return the first tile URL that contains '.pbf' (even with query strings).
    Supports relative tile URLs by resolving against the tiles.json base.
    """
    if url.endswith(".json") or "tiles.json" in url:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        meta = r.json()
        base = url.rsplit("/", 1)[0] + "/"
        tiles = meta.get("tiles", [])
        for t in tiles:
            if isinstance(t, str) and ".pbf" in t and "{z}" in t and "{x}" in t and "{y}" in t:
                return t if t.startswith("http") else urljoin(base, t)
            if isinstance(t, dict):
                candidate = t.get("url") or t.get("href") or ""
                if ".pbf" in candidate and "{z}" in candidate and "{x}" in candidate and "{y}" in candidate:
                    return candidate if candidate.startswith("http") else urljoin(base, candidate)
        raise RuntimeError("Could not find a .pbf template in tiles.json")
    else:
        return url

def http_or_cache(url: str) -> bytes:
    key = hashlib.sha1(url.encode("utf-8")).hexdigest()[:40]
    cache_path = CACHE_DIR / f"{key}.pbf"
    if cache_path.exists():
        return cache_path.read_bytes()
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    cache_path.write_bytes(r.content)
    return r.content

def read_tile_http(pbf_template: str, z: int, x: int, y: int) -> bytes:
    return http_or_cache(pbf_template.format(z=z, x=x, y=y))

def read_tile_from_mbtiles(mbtiles_path: str, z: int, x: int, y: int) -> Optional[bytes]:
    import sqlite3
    if not mbtiles_path:
        return None
    if not os.path.exists(mbtiles_path):
        raise FileNotFoundError(f"MBTiles not found: {mbtiles_path}")
    conn = sqlite3.connect(mbtiles_path)
    try:
        cur = conn.cursor()
        tms_y = (2 ** z - 1) - y
        cur.execute("SELECT tile_data FROM tiles WHERE zoom_level=? AND tile_column=? AND tile_row=?", (z, x, tms_y))
        row = cur.fetchone()
        return row[0] if row and row[0] else None
    finally:
        conn.close()

def load_mvt_polygons_for_tile(raw_bytes: bytes, z: int, x: int, y: int, layer_allowlist: List[str]) -> List[shapely.geometry.base.BaseGeometry]:
    decoded = mapbox_vector_tile.decode(raw_bytes)

    if not layer_allowlist:
        print(f"[INFO] Layers in tile z{z}/{x}/{y}: {list(decoded.keys())}")
        candidate_layers = list(decoded.keys())
    else:
        candidate_layers = [lyr for lyr in decoded.keys() if lyr in layer_allowlist]

    polys: List[shapely.geometry.base.BaseGeometry] = []
    for lyr in candidate_layers:
        features = decoded[lyr]["features"]
        extent = decoded[lyr].get("extent", 4096)
        for feat in features:
            props = feat.get("properties", {})
            if PROPERTY_FILTER:
                ok = True
                for k, allowed_vals in PROPERTY_FILTER.items():
                    if str(props.get(k)) not in set(allowed_vals):
                        ok = False
                        break
                if not ok:
                    continue
            geom = mvt_geom_to_webm(z, x, y, feat["geometry"], extent=extent)
            if geom is None:
                continue
            if isinstance(geom, (Polygon, MultiPolygon)):
                polys.append(geom)
    return polys

def max_continuous_overlap_m(line_webm, danger_union) -> float:
    if danger_union.is_empty:
        return 0.0
    inter = line_webm.intersection(danger_union)
    if inter.is_empty:
        return 0.0
    lengths = []
    if isinstance(inter, LineString):
        lengths.append(inter.length)
    elif isinstance(inter, MultiLineString):
        for seg in inter.geoms:
            lengths.append(seg.length)
    elif isinstance(inter, GeometryCollection):
        for g in inter.geoms:
            if isinstance(g, (LineString, MultiLineString)):
                if isinstance(g, LineString):
                    lengths.append(g.length)
                else:
                    for seg in g.geoms:
                        lengths.append(seg.length)
    return max(lengths) if lengths else 0.0

def classify_route(max_cont_m: float) -> str:
    if max_cont_m >= THRESH_RED_MIN_CONTINUOUS:
        return "red"
    elif max_cont_m < THRESH_GREEN_MAX_CONTINUOUS:
        return "green"
    else:
        return "orange"

def covering_tiles_for_geom(geom_wgs84, z: int):
    minx, miny, maxx, maxy = geom_wgs84.bounds
    for t in mercantile.tiles(minx, miny, maxx, maxy, [z]):
        yield (t.z, t.x, t.y)

def process_level(level_name: str, tileset_url: str, gdf: gpd.GeoDataFrame):
    print(f"\n=== Processing level: {level_name} ===")
    pbf_template = resolve_pbf_template(tileset_url)

    col_len = f"max_continuous_{level_name}_m"
    col_cls = f"ass_class_{level_name}"
    if col_len not in gdf.columns: gdf[col_len] = 0.0
    if col_cls not in gdf.columns: gdf[col_cls] = None

    for idx, row in gdf.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            gdf.at[idx, col_len] = 0.0
            gdf.at[idx, col_cls] = classify_route(0.0)
            continue

        tiles = sorted(set(covering_tiles_for_geom(geom, Z)))
        if not tiles:
            gdf.at[idx, col_len] = 0.0
            gdf.at[idx, col_cls] = classify_route(0.0)
            continue

        danger_polys_webm = []
        for (tz, tx, ty) in tiles:
            raw = None
            mb = READ_FROM_MBTILES.get(level_name, "")
            if mb:
                raw = read_tile_from_mbtiles(mb, tz, tx, ty)
            if raw is None:
                raw = read_tile_http(pbf_template, tz, tx, ty)
            polys = load_mvt_polygons_for_tile(raw, tz, tx, ty, DANGER_LAYER_ALLOWLIST)
            danger_polys_webm.extend(polys)

        if not danger_polys_webm:
            gdf.at[idx, col_len] = 0.0
            gdf.at[idx, col_cls] = classify_route(0.0)
            continue

        danger_union = unary_union(danger_polys_webm)
        line_webm = shapely.ops.transform(PROJECT_TO_WEBM, geom)
        max_cont = max_continuous_overlap_m(line_webm, danger_union)

        gdf.at[idx, col_len] = float(max_cont)
        gdf.at[idx, col_cls] = classify_route(float(max_cont))
        print(f"[{idx}] {level_name}: max_cont={max_cont:.2f} m -> {gdf.at[idx,col_cls]}")

def main():
    if not os.path.exists(SKI_ROUTES_PATH):
        print(f"ERROR: Can't find ski routes GeoJSON at {SKI_ROUTES_PATH}")
        sys.exit(2)

    gdf = gpd.read_file(SKI_ROUTES_PATH)
    if gdf.crs is None:
        print("[WARN] Input has no CRS; assuming EPSG:4326 (WGS84 lon/lat).")
        gdf = gdf.set_crs(4326)
    elif gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs(4326)

    level_order = ["low","moderate","scary_moderate","considerable","high"]
    for lvl in level_order:
        tileset_url = DEFAULT_DANGER_TILESETS.get(lvl)
        if not tileset_url:
            print(f"[WARN] No tileset URL for level '{lvl}', skipping.")
            continue
        process_level(lvl, tileset_url, gdf)

    gdf.to_file(OUTPUT_GEOJSON, driver="GeoJSON")
    print(f"\nDone. Wrote {OUTPUT_GEOJSON}")
    print("Per-level columns added:")
    for lvl in level_order:
        print(f"  - max_continuous_{lvl}_m")
        print(f"  - ass_class_{lvl} (red/orange/green)")

    if not DANGER_LAYER_ALLOWLIST:
        print("\nTIP: Re-run with env var DANGER_LAYERS='layer1,layer2' once you know the layer name(s).")

if __name__ == "__main__":
    main()
