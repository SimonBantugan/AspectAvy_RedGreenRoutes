# AspectAvy_RedGreenRoutes
Code for red / green qualification for routes based on Adaptive Slope Shading 
# 🏔️ Ski Route Avalanche Danger Classifier (MapTiler Vector Tiles)

This tool measures how much of each ski line overlaps **avalanche danger polygons** from MapTiler vector tiles.  
It works for **five danger levels** — *low*, *moderate*, *scary moderate*, *considerable*, and *high* —  
and classifies each route segment based on **continuous avalanche slope stability (ASS) danger length**.

---

## ✨ Key features
- Streams only the tiles that overlap each ski route — no massive downloads.
- Supports both **online MapTiler vector tiles** and **local `.mbtiles`**.
- Computes the **maximum continuous overlap length** (in meters) of each route with avalanche polygons.
- Adds per-level attributes:

  | Column | Description |
  |--------|--------------|
  | `max_continuous_<level>_m` | Longest continuous ASS-danger overlap (m) |
  | `ass_class_<level>` | `red` ≥ 40 m, `green` < 12 m, `orange` otherwise |

- Fully configurable thresholds, cache directory, and zoom level.

---

## 📦 Requirements
Python 3.9 +   
Dependencies:
```bash
pip install geopandas shapely pyproj mercantile requests mapbox-vector-tile
