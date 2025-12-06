import pandas as pd
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import os
import sys

colors = ['#8EC8E2', '#E47B81', '#60B67B', '#FF8C42', '#F2C94C', '#B05FA0']

sample_list = ["ERR10513574","ERR10963128",
                 "ERR11026612","ERR11026613","ERR11026615","ERR11026616","ERR11026617","ERR11026625","ERR11026630","ERR11026631","ERR11026635","ERR11026636","ERR11026637","ERR11026638",
                 "ERR11030164","ERR11030165","ERR11030167","ERR11030168","ERR11030169","ERR11030177","ERR11030182","ERR11030183","ERR11030187","ERR11030188","ERR11030189","ERR11030190",
                 "SRR22450503","SRR22450504","SRR22450505","SRR22450506","SRR22450507","SRR22450508","SRR22450509","SRR22450510","SRR22450511",
                 "SRR22450515","SRR22450516","SRR22450517","SRR22450518","SRR22450519","SRR22450520","SRR22450521","SRR22450522","SRR22450523",
                 "SRR31266107","SRR31266108","SRR31266109","SRR31266110","SRR31266111","SRR31266112","SRR31266113",
                 "SRR31266114","SRR31266115","SRR31266116","SRR31266117","SRR31266118","SRR31266119","SRR31266120",
                 "SRR19536726","SRR19536727","SRR19536728","SRR19536729",
                 "SRR24877167","SRR24877168"]

project_list = ["PRJEB56841","PRJEB56841",
                  "PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728",
                  "PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728",
                  "PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618",
                  "PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618",
                  "PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318",
                  "PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318",
                  "PRJNA845087","PRJNA845087","PRJNA845087","PRJNA845087",
                  "PRJNA981509","PRJNA981509"]

for i,sample in enumerate(sample_list):

 project = project_list[i]

 for caller in ["clair3rna","lofreq","bcftools"]:

    df = pd.read_csv("data/"+sample+"_"+caller+".txt",sep='\t')
    #df = df.rename(columns={"V5": "alt"})
    #print(df.size)
    if project == "PRJEB60728":
     df = df[df.vaf >= 0.1]
    else:
     df = df[df.vaf >= 0.01]

    #print(df.size)
    #sys.exit(1)
    dfgrp = df.groupby(["motif3","alt"])["altcnt"].sum().reset_index()
    #dfgrp = df.groupby(["motif3","alt"]).size().reset_index(name='Count')

    A = ['A'.join(i) for i in itertools.product('ACGT', repeat=2)]
    C = ['C'.join(i) for i in itertools.product('ACGT', repeat=2)]
    G = ['G'.join(i) for i in itertools.product('ACGT', repeat=2)]
    T = ['T'.join(i) for i in itertools.product('ACGT', repeat=2)]

    mt = pd.DataFrame(list(itertools.product(A+C+G+T,'ACGT')), columns=["motif3","alt"])
    mt["ref"] = mt.motif3.str[1]
    mt = mt[mt.ref != mt.alt]

    mt1 = mt[(mt.motif3.str[1] == 'C') | (mt.motif3.str[1] == 'T')]
    mt1 = mt1.sort_values(by=["ref","alt","motif3"]).reset_index()
    mt2 = mt[(mt.motif3.str[1] == 'A') | (mt.motif3.str[1] == 'G')]
    mt2 = mt2.sort_values(by=["ref","alt","motif3"]).reset_index()

    plotdf1 = pd.merge(mt1, dfgrp, on=["motif3","alt"], how="left")
    plotdf1["altcnt"].fillna(0, inplace=True)
    #plotdf1["Count"].fillna(0, inplace=True)

    plotdf2 = pd.merge(mt2, dfgrp, on=["motif3","alt"], how="left")
    plotdf2["altcnt"].fillna(0, inplace=True)
    #plotdf2["Count"].fillna(0, inplace=True)

    fig = plt.figure(figsize=(16, 7))

    ax1 = plt.subplot2grid(shape=(2, 1), loc=(0, 0))
    ax2 = plt.subplot2grid(shape=(2, 1), loc=(1, 0))

    custom_palette = []
    for x in colors:
        custom_palette.extend([x]*16)
    ax1 = sns.barplot(x=plotdf1.index, y=plotdf1.altcnt, palette=custom_palette, ax=ax1)
    #ax1 = sns.barplot(x=plotdf1.index, y=plotdf1.Count, palette=custom_palette, ax=ax1)
    ax2 = sns.barplot(x=plotdf2.index, y=plotdf2.altcnt, palette=custom_palette, ax=ax2)
    #ax2 = sns.barplot(x=plotdf2.index, y=plotdf2.Count, palette=custom_palette, ax=ax2)
    ax1.set_xticklabels(plotdf1.motif3, rotation=90, horizontalalignment='center')
    ax2.set_xticklabels(plotdf2.motif3, rotation=90, horizontalalignment='center')

    ymax1 = ax1.get_ylim()[1]+ax1.get_ylim()[1]/10
    ax1.set_ylim(0, ymax1)
    ax1.set_xlabel("Motif")
    ax1.set_ylabel("Number of substitutions")
    ymax2 = ax2.get_ylim()[1]+ax2.get_ylim()[1]/10
    ax2.set_ylim(0, ymax2)
    ax2.set_xlabel("Motif")
    ax2.set_ylabel("Number of substitutions")

    sbs1 = ["C>A","C>G","C>T","T>A","T>C","T>G"]
    sbs2 = ["A>C","A>G","A>T","G>A","G>C","G>T"]
    for i in range(len(colors)):
        ax1.add_patch(Rectangle((i*16, ymax1-ymax1/10), 16, ymax1/10, facecolor = colors[i]))
        ax1.text(i*16+8, ymax1-ymax1/20, sbs1[i], horizontalalignment='center', verticalalignment='center', size='x-large', color='white', weight='semibold')
        ax2.add_patch(Rectangle((i*16, ymax2-ymax2/10), 16, ymax2/10, facecolor = colors[i]))
        ax2.text(i*16+8, ymax2-ymax2/20, sbs2[i], horizontalalignment='center', verticalalignment='center', size='x-large', color='white', weight='semibold')

    plt.tight_layout()
    plt.close(fig)

    save_path = "PICS/fig1/c/"+project+"/fig1c_motif3_"+caller+"_"+project+"_"+sample+".png"
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    fig.savefig(save_path)
