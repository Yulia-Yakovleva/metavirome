import folium
import pandas as pd
from folium.plugins import MarkerCluster

data = pd.read_csv(r"coordinates.csv", na_values=["None"])
data = data.dropna()

m = folium.Map()
cluster = MarkerCluster().add_to(m)

for d in data.itertuples():
    folium.Marker([d.latitude, d.longitude], popup=d.id2, tooltip=d.id).add_to(cluster)

m.save("map.html")
