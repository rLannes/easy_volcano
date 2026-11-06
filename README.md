# easy_volcano
Romain JSB Lannes, PhD

Fully customisable publication-ready volcano plots in one line.


easy_volcano are a plots used ubiquitously in RNA-seq analyses, yet in python there are few alternatives than to write custom code each time. This library aims is to allow make highly customizable volcano-plot in python in one line. 

# Installation
pip install easy_volcano

# Quick start
simpliest example:

```
import pandas as pd
import easyvolcano

example_file= "" # in this folder
df = pd.read_csv(example_file, sep='\t')
(fig, ax) = easyvolcano.plot_volcano(df)
plt.plot()

```

# Almost endless personalization 

easy_volcano allow very deep personalization by using a (selector, override). 
the easiest way to understand it to look at the example:

Where:
    - Selector allow to select genes to apply specific visual personalization. It can be a gene, a list of gene, a mask, a lamdba function. 
    - Override must be a dict that accept:
        - **kwargs argument for scatter plot
        - label: a string, will add it to the legend.
        - label_gene: bolean or dict, if dict pass **kwargs to modify text

# Example

```
(fig, ax) = easyvolcano.plot_volcano(df, )
plt.plot()
```