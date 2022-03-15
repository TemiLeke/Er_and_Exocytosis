import random
import numpy as np
import math
import matplotlib.pyplot as plt
import time
import pandas as pd
from sklearn.linear_model import LinearRegression


x = []
y = []

for line in open("ReleaseRate.csv"):
    x.append(line.split()[0])
    y.append(line.split()[0])

x = [float(i) for i in x]
y = [float(j) for j in y]

my_df = pd.DataFrame({'x': x , 'y': y})
my_df.loc[:, 'log_x'] = map(lambda x: math.log(x), my_df['x'])
my_df.loc[:, 'log_y'] = map(lambda x: math.log(x), my_df['y'])
model = LinearRegression()
model.fit(my_df['log_x'].values.reshape(-1, 1), my_df['log_y'].values)
print(model.coef_)
print(model.intercept_)

X = 0.25
log_y = model.intercept_ + model.coef_ * math.log(X)
y = math.exp(log_y)
print(y)
