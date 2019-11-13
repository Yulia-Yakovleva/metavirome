import folium
import pandas as pd

data = pd.read_csv(r"coordinates.csv", na_values=["None"])
data = data.dropna()

m = folium.Map()

for d in data.itertuples():
    folium.Marker([d.latitude, d.longitude], popup=d.id, tooltip=d.id2).add_to(m)

m.save("map.html")
