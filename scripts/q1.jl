using JackknifeAnalysis
using Statistics

v1 = readdata("data/v1")
v2 = readdata("data/v2")
v3 = readdata("data/v3")
v4 = readdata("data/v4")
v5 = readdata("data/v5")

map(mean, (v1, v2, v3, v4, v5))
