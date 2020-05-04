from fyodor import pwv, cerro_paranal
import matplotlib.pyplot as plt

date, water_vapor = pwv(cerro_paranal, Ra=200, Dec=50)

print(date, water_vapor)
