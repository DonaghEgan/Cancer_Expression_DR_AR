# Title: plot meta data figure
# Description:  read in clincal annotations and plot. Code is adopted and modified
# from Liu et al 
# Author: Donagh
# Date: 2025-01-28
# Version: 1.1

import pyreadr
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
import pandas as pd
import matplotlib.patches as mpatches
import numpy as np
import brewer2mpl

brewer2mpl.print_maps()
COLOR_LIST = brewer2mpl.get_map('Set3', 'Qualitative', 12).mpl_colors

GRAY = (0.95, 0.95, 0.95, 1)
LIGHT_BLUE = (0.65098041296005249, 0.80784314870834351, 0.89019608497619629, 1.0)
BLUE = (0.12572087695201239, 0.47323337360924367, 0.707327968232772, 1.0)
LIGHT_GREEN = (0.68899655751153521, 0.8681737867056154, 0.54376011946622071, 1.0)
GREEN = (0.21171857311445125, 0.63326415104024547, 0.1812226118410335, 1.0)
PINK = (0.98320646005518297, 0.5980161709820524, 0.59423301088459368, 1.0)
RED = (0.89059593116535862, 0.10449827132271793, 0.11108035462744099, 1.0)
LIGHT_ORANGE = (0.99175701702342312, 0.74648213716698619, 0.43401768935077328, 1.0)
ORANGE = (0.99990772780250103, 0.50099192647372981, 0.0051211073118098693, 1.0)
LIGHT_PURPLE = (0.78329874347238004, 0.68724338552531095, 0.8336793640080622, 1.0)
PURPLE = (0.42485198495434734, 0.2511495584950722, 0.60386007743723258, 1.0)
YELLOW = (0.99760092286502611, 0.99489427150464516, 0.5965244373854468, 1.0)
LIGHT_BROWN = (0.9, 0.6, 0.5, 1.0)
BROWN = (0.69411766529083252, 0.3490196168422699, 0.15686275064945221, 1.0)

YELLOW = (255/256.,255/256.,51/256.,1.0)
BLUE = (55/256.,126/256.,184/256.,1.0)
PURPLE = (152/256.,78/256.,163/256.,1.0)
RED = (228/256.,26/256.,28/256.,1.0)
GREEN = (77/256.,175/256.,74/256.,1.0)
BROWN = (166/256.,86/256.,40/256.,1.0)
ORANGE = (255/256.,127/256.,0,1.0)
GRAY = (0.95, 0.95, 0.95, 1.0)

RESPONDER_COL = "forestgreen"
NRESPONDER_COL = "goldenrod"

NF_COL = "#4dbbd5"
AR_COL = "#e64b35"
PR_COL = "#00a087"
NA_COL = GRAY

SKIN_COL = BLUE
OTHER_COL = YELLOW
LN_COL = PINK
LUNG_COL = LIGHT_GREEN
LIV_COL = LIGHT_PURPLE

PD1_COL = "cornflowerblue"
COMBO_COL = LIGHT_GREEN
IPI_COL = "goldenrod"

R_patch = mpatches.Patch(color = RESPONDER_COL, label = "Response")
NR_patch = mpatches.Patch(color = NRESPONDER_COL, label = "NonResponse")

NF_patch = mpatches.Patch(color = NF_COL, label = "Nofailure")
AR_patch = mpatches.Patch(color = AR_COL, label = "Acq_Resistance")
PR_patch = mpatches.Patch(color = PR_COL, label = "Primary_Resistance") 
NA_patch = mpatches.Patch(color = NA_COL, label = "Not Assigned") 

SKIN_patch = mpatches.Patch(color = SKIN_COL, label = "Skin")
OTHER_patch = mpatches.Patch(color = OTHER_COL, label = "Other")
LN_patch = mpatches.Patch(color = LN_COL, label = "LN")
LIVER_patch = mpatches.Patch(color = LIV_COL, label = "Liver")
LUNG_patch = mpatches.Patch(color = LUNG_COL, label = "Lung")

PD1_patch = mpatches.Patch(color = PD1_COL, label = "anti-PD1")
COMBO_patch = mpatches.Patch(color = COMBO_COL, label = "Combination")
IPI_patch = mpatches.Patch(color = IPI_COL, label = "Ipilimumab")

Alive_patch = mpatches.Patch(color = "#fde725", label = "Alive")
Dead_patch = mpatches.Patch(color = "#440154", label = "Dead")

meta_path =  "data/meta/meta_data.Rds" # Replace with the actual file path
counts_path = "available from HMF" # Replace with the actual file path

# Read the RDS file
result = pyreadr.read_r(meta_path)
counts = pyreadr.read_r(counts_path)

# Access the data frame from the result
counts = counts[None]
data_frame = result[None]  # None is used to access the data frame
data_frame = data_frame[data_frame.index.isin(counts.columns)]

# reassign sites
data_frame.loc[data_frame['BiopsySiteClean'] == 'IT?', 'BiopsySiteClean'] = 'Other'
data_frame.loc[data_frame['BiopsySiteClean'] == 'Bone', 'BiopsySiteClean'] = 'Other'
data_frame.loc[data_frame['BiopsySiteClean'] == 'GI', 'BiopsySiteClean'] = 'Other'

data_frame.loc[data_frame['BiopsySiteClean'] == 'Skin', 'BiopsySiteClean'] = 'Skin'
data_frame.loc[data_frame['BiopsySiteClean'] == 'SoftTissue', 'BiopsySiteClean'] = 'Skin'
data_frame.loc[data_frame['BiopsySiteClean'] == 'Subcut', 'BiopsySiteClean'] = 'Skin'

##
order = data_frame.sort_values(by=['group', 'Duration'], ascending=False).index
data_frame.loc[order, 'FirstIOResNoR'] = data_frame.loc[order, 'FirstIOResNoR'].map({'NonResponse': 0, 'Response': 1})
data_frame.loc[order, 'group'] = data_frame.loc[order,'group'].map({'Primary_Resistance': 0, 'Acq_Resistance': 1, 'Nofailure': 2, 'NA': 3})
data_frame.loc[order, 'BiopsySiteClean'] = data_frame.loc[order,'BiopsySiteClean'].map({'Other': 0, 'Skin': 1, 'Lung': 2, 'LN':3, 'Liver':4})
data_frame.loc[order, 'treatment'] = data_frame.loc[order,'treatment'].map({'anti-PD1': 0, 'Ipilimumab': 1, 'Combination': 2})
data_frame.loc[order, 'VitalStatus'] = data_frame.loc[order,'VitalStatus'].map({'Alive': 0, 'Dead': 1})
data_frame = data_frame.T

#data_frame['response_bin'] = data_frame['response_bin'].replace({'Non_Responder': 0, 'Responder': 1})
#data_frame['group'] = data_frame['group'].replace({'Primary_Resistance': 0, 'Acq_Resistance': 1, 'Nofailure': 2})
#data_frame = data_frame.T
#order = data_frame.sort_values(by=['tumorPurity', 'Duration'], ascending=False).index

fig = fig = plt.figure(figsize=(11, 6), constrained_layout=True)
sizeAdj = 0.6


nrows = 7
height_ratios = [4,1,1,1,1,1,1]
ncols = 3
    
gs = gridspec.GridSpec(nrows,ncols, height_ratios = height_ratios, width_ratios = [5,1,1], 
                           wspace=0.1 * sizeAdj, hspace = 0.05* sizeAdj)
    
coMutPlot = plt.subplot(gs[6,0])

def PlotStackedBarPlot(barPlot, df_nonnorm, sigOrder = None, xlabel = False, 
                       normalize = False, legendSize = 50):
    numSamples = len(df_nonnorm.columns)
    numClasses = len(df_nonnorm)

#    plt.figure(figsize=(30,10))

    if normalize: 
        # normalize data
        df = df_nonnorm/df_nonnorm.apply(sum, axis = 0)
    else:
        df = df_nonnorm

    # set order
    if sigOrder is None: # set signature order to default order
        sigOrder = df_nonnorm.index    
    df = df.loc[sigOrder]
    
    # generate the bottoms of each of the stacked values
    bottom = [np.array([0] * numSamples)]
    for i in range(0, numClasses):
        bottom.append(bottom[i] + np.array(df.iloc[i,:]))
    
    ind = np.arange(numSamples)
    width = 0.5
    bargraphs = []
    
    for i in range(0, numClasses):
        bargraphs.append(plt.bar(ind, df.iloc[i,:], width, bottom = bottom[i]))

    # label x axis
    if xlabel:
        plt.setp(barPlot.get_xticklabels(), visible=True)
        plt.xticks(ind, df.columns, rotation = 'vertical')
#        barPlot.get_xaxis().set_ticks(ind + width/2., df.columns, rotation = 'vertical')
    else:
        plt.setp(barPlot.get_xticklabels(), visible=False)

    # Limit x axis
    plt.xlim(-0.5, numSamples-0.5)
    
#        barPlot.get_xaxis().set_ticks(ind + width/2., [])
    
    # label y axis
    barPlot.get_yaxis().set_label("Mutations")
    
    
    # Legend
    if legendSize is not None:
        barPlot.legend([x[0] for x in bargraphs[::-1]], sigOrder[::-1], prop = {'size':legendSize}, 
                        bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

# example
testDF = pd.DataFrame({'sample4': [87,37,0,5],
                       'sample1': [35,50,10,5],
                       'sample2': [15,3,50,5],
                       'sample3': [30,5,20,1]},
                      index=['sig1','sig2','sig3','sig4'],
                      columns=['sample4','sample1','sample2','sample3'])
                       
fig = plt.figure(figsize = (5,5))
                       
def CreateCoMutFig(currentPlot, df, colorList = COLOR_LIST, cmap = None, italic=False):
    numGenes = len(df)
    numSamples = len(df.columns)

    # remove initial frame
    for side in ['top', 'bottom', 'right', 'left']: 
        currentPlot.spines[side].set_visible(False)
    
    # set colors
    if cmap is None:
        cmap = mpl.colors.ListedColormap(colorList)
        vmax = len(colorList) - 1
    else:
        vmax = 1
    
    # display data
    cax = currentPlot.imshow(df, interpolation = 'none', aspect = 'auto',
                    cmap = cmap, vmin = 0, vmax = vmax)
    
    # Remove tick marks
    plt.tick_params(axis = 'both', which = 'both', 
                top = False, bottom = False, right = False, left = False)
    
    # Resize tick marks
    plt.xticks(range(0, numSamples, 1), rotation = 90)
    plt.yticks(range(0, numGenes,1), fontsize = 13)
    
    # Replace y ticks with gene names
    if italic:
        currentPlot.set_yticklabels(df.index, style = 'italic')
    else:
        currentPlot.set_yticklabels(df.index)

    currentPlot.set_xticklabels(df.columns)
 
    
    # Create white lines across the grid
    for y in range(0, numGenes, 1):  
        plt.plot(range(-1, numSamples+2), [y+0.5] * (len(range(-1, numSamples))+2), "-", 
                 lw=2, color="white", alpha=1)
    for x in range(0,numSamples+1,1):
        plt.plot([x + 0.5] * len(range(-1,numGenes+1)), range(-1,numGenes+1), "-",
                 lw=2, color="white", alpha=1)    
        
    # limit the x and y-axis to just the range
    plt.xlim(-0.5, numSamples-0.5)
    plt.ylim(-0.5, numGenes-0.5)
    
    # remove grid lines
    currentPlot.grid(False)
    return (cax)

#generate group plot
groupPlot = plt.subplot(gs[1,0], sharex = coMutPlot)
CreateCoMutFig(groupPlot, data_frame.loc[['group'],order].astype(int), [PR_COL, AR_COL, NF_COL, NA_COL])
plt.setp(groupPlot.get_xticklabels(), visible = False) # suppress x labels
plt.setp(groupPlot.get_yticklabels(), size=20 * sizeAdj) 

#generate response plot
responsePlot = plt.subplot(gs[2,0], sharex = coMutPlot)
CreateCoMutFig(responsePlot, data_frame.loc[['FirstIOResNoR'],order].astype(int), [NRESPONDER_COL, RESPONDER_COL])
plt.setp(responsePlot.get_xticklabels(), visible = False) # suppress x labels
plt.setp(responsePlot.get_yticklabels(), size=20 * sizeAdj) 

#generate purity plot
purityPlot = plt.subplot(gs[3,0], sharex = coMutPlot)
CreateCoMutFig(purityPlot, data_frame.loc[['tumorPurity'],order].astype(float),  cmap=plt.get_cmap("Purples"))
plt.setp(purityPlot.get_xticklabels(), visible = False) # suppress x labels
plt.setp(purityPlot.get_yticklabels(), size=20 * sizeAdj) 

#generate biopsy plot
biopsyPlot = plt.subplot(gs[4,0], sharex = coMutPlot)
CreateCoMutFig(biopsyPlot, data_frame.loc[['BiopsySiteClean'],order].astype(int),  [LN_COL, SKIN_COL, LUNG_COL, LIV_COL, ])
plt.setp(biopsyPlot.get_xticklabels(), visible = False) # suppress x labels
plt.setp(biopsyPlot.get_yticklabels(), size=20 * sizeAdj) 

#generate treatment plot
TrxPlot = plt.subplot(gs[5,0], sharex = coMutPlot)
CreateCoMutFig(TrxPlot, data_frame.loc[['treatment'],order].astype(int),  [PD1_COL, COMBO_COL, IPI_COL])
plt.setp(TrxPlot.get_xticklabels(), visible = False) # suppress x labels
plt.setp(TrxPlot.get_yticklabels(), size=20 * sizeAdj) 

#generate os plot
OSPlot = plt.subplot(gs[6,0], sharex = coMutPlot)
CreateCoMutFig(OSPlot, data_frame.loc[['VitalStatus'],order].astype(int),  ["#fde725", "#440154"])
plt.setp(OSPlot.get_xticklabels(), visible = False) # suppress x labels
plt.setp(OSPlot.get_yticklabels(), size=20 * sizeAdj) 

#generate follow up on top
mutLoadBarPlot = plt.subplot(gs[0,0], sharex = coMutPlot)
PlotStackedBarPlot(mutLoadBarPlot, data_frame.loc[['Duration'], order],
                   normalize=False, legendSize = 0) 
mutLoadBarPlot.set_ylabel("Follow Up (Days)", fontsize=20 * sizeAdj)
plt.axhline(y=182.5, color='red', linestyle='--', label='Min')
plt.setp(mutLoadBarPlot.get_yticklabels(), size=15 * sizeAdj) 
plt.ticklabel_format(style='plain', axis='y')


###################
# LEGENDS
###################

#Group legend
curPlot = plt.subplot(gs[0,1])
legend = curPlot.legend(handles = [NF_patch, AR_patch, PR_patch, NA_patch], title="Group",
              prop = {'size':15* sizeAdj}, bbox_to_anchor=(0,1), loc="upper left", borderaxespad=0.)
plt.setp(legend.get_title(),fontsize=15* sizeAdj)
plt.axis("off")
plt.setp(coMutPlot.get_xticklabels(), visible = False)

#Response legend
curPlot = plt.subplot(gs[0,2])
legend = curPlot.legend(handles = [R_patch, NR_patch], title="Response",
              prop = {'size':15* sizeAdj}, bbox_to_anchor=(0,1), loc="upper left", borderaxespad=0.)
plt.setp(legend.get_title(),fontsize=15* sizeAdj)
plt.axis("off")
plt.setp(coMutPlot.get_xticklabels(), visible = False)

#Biopsy legend
curPlot = plt.subplot(gs[1,1])
legend = curPlot.legend(handles = [LN_patch, SKIN_patch, OTHER_patch, LIVER_patch, LUNG_patch], title="Biopsy Site",
              prop = {'size':15* sizeAdj}, bbox_to_anchor=(0,1), loc="upper left", borderaxespad=0.)
plt.setp(legend.get_title(),fontsize=15* sizeAdj)
plt.axis("off")
plt.setp(coMutPlot.get_xticklabels(), visible = False)

#trx legend
curPlot = plt.subplot(gs[1,2])
legend = curPlot.legend(handles = [PD1_patch, COMBO_patch, IPI_patch], title="Treatment",
              prop = {'size':15* sizeAdj}, bbox_to_anchor=(0,1), loc="upper left", borderaxespad=0.)
plt.setp(legend.get_title(),fontsize=15* sizeAdj)
plt.axis("off")
plt.setp(coMutPlot.get_xticklabels(), visible = False)

#os legend
curPlot = plt.subplot(gs[2,2])
legend = curPlot.legend(handles = [Alive_patch, Dead_patch], title="OS",
              prop = {'size':15* sizeAdj}, bbox_to_anchor=(0,1), loc="upper left", borderaxespad=0.)
plt.setp(legend.get_title(),fontsize=15* sizeAdj)
plt.axis("off")
plt.setp(coMutPlot.get_xticklabels(), visible = False)

plt.show()
