# multiomics_correlation_plot
Merges and correlates multiomics data and creates visualizations with significantly correlated features. 

# What data can I use?
Anything that is a continuous variable. In the example here, I compare untargeted metabolomics, 16S microbiome sequencing, and human milk oligosaccharide (HMO) data. However, before you compare datasets, it's important to consider what type of data you're comparing. Data sparseness, normalization, and missing data imputation can all affect the outcomes of multi-omics comparisons. Here, data was processed using rclr and missing data imputation before running this script (for more information, see https://github.com/Sydney-Thomas/iimn_clean_normalize). At the very least, I suggest log centering data before running correlations. 

# What is being compared?
You choose! In the example data provided, I'm interested in looking at correlations between 3 features (2 microbes and 1 HMO) and the entire metabolome. To do so, I've created a dataframe with those three "features of interest" in the "From" column, and all the metabolites in the "To" column. These comparisons can be whatever you want, just make sure the comparisons you're interested in are included in those two columns.   

On the statistical side, Spearman correlations are used since they don't assume normality. p values are then adjusted using Benjamini-Hochberg and only features that are significantly correlated to all three "features of interest" are kept. All of these parameters can be adjusted.

# What are the outputs?
Two correlation plots are created. The first is a circle correlation plot. This type of plot is useful when you want an overview of the data without any prior knowledge of what factors are important.
![Circle](https://github.com/Sydney-Thomas/multiomics_correlation_plot/blob/1b685ccaf1a4236cbf4caffebaee2ed1b4e7a4f4/Circle_Correlation.png)

The next is a column correlation plot that highlights specific features. This type of plot is useful if you already have a handful of features that you're interested in.
![Column](https://github.com/Sydney-Thomas/multiomics_correlation_plot/blob/1b685ccaf1a4236cbf4caffebaee2ed1b4e7a4f4/Column_Correlation.png)
