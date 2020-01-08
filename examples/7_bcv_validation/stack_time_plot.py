
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set(style="ticks", color_codes=True)
df = pd.read_csv("cav_time.csv")

sns.catplot(kind='bar', x='Method', y='Time', col='Step', data=df, palette='Reds', sharey=False, ci=None, col_wrap=2)

plt.savefig("time.png")