import pandas as pd

data = pd.DataFrame()
for file in ["coordinates_marine.csv", "coordinates_sea.csv", "coordinates_soil.csv"]:
    data = pd.concat([data, pd.read_csv(file)])

print(len(data))
data = data.drop_duplicates()
print(len(data))
data.to_csv("coordinates.csv")
