import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


def scatter(df):

    sns.set(rc = {'figure.figsize': (20,10)})
    sns.set(font_scale=5 )

    g = sns.FacetGrid(df, col= "target", col_wrap=3, height=15, aspect=1, legend_out=True, margin_titles=True,
                      gridspec_kws={"wspace":0.18})

    g = g.map(sns.scatterplot, 'x', 'y', 'ligand', 'score', palette='Paired', alpha=0.9, s=2200, edgecolor='b')


    g.map(sns.lineplot, x=[0,100], y =[0, 100], color='grey')

    plt.savefig("scatter4.png")
    plt.close()


def main():
    sns.set_style("dark")
    sns.set_context("paper")
    df = pd.read_csv('hot_vs_bcv.csv')
    print(df)

    scatter(df)



if __name__ == '__main__':
    main()
