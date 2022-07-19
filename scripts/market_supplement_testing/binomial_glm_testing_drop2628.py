import geopandas as gpd
import pandas as pd 
import numpy as np
from geopy import distance
import statsmodels.formula.api as smf
import statsmodels.api as sm1

gdf = gpd.read_file("../../maps/geojson/huanan-market-internal.geojson")

pos_sites =gdf[gdf['group']=='Env-Pos']
neg_stalls =gdf[gdf['group']=='Env-Neg']

## get location of environmental negatives 
series_neg = gpd.GeoSeries(neg_stalls.loc[:,'geometry'])
df_neg = pd.DataFrame()
df_neg['x'] = series_neg.x
df_neg['y'] = series_neg.y
df_neg = df_neg.reset_index(drop=True)

## get location of environmental positives 
pos_sites['title'] = pos_sites['title'].apply(lambda x: x.split(' ')[0]) #strips sample type info
pos_sites = pos_sites.set_index('title')
series = gpd.GeoSeries(pos_sites.loc[:,'geometry'])
df_pos = pd.DataFrame()
df_pos['x'] = series.x
df_pos['y'] = series.y

# Stalls with multiple environmental positives, using more specific placement from CCDC map.
# Locations are centroids of each business, for positivity analysis. 
same_stall=[['A14','A15','A90','A88','A87'],['Q61','Q64','Q69','Q70','Q68'],['A2','A18','A20'],['F98','F100'],['G93','A61']]
same_stall_center = [[114.25658,30.61937],[114.256469,30.61947],[114.256724,30.619613],[114.25704,30.61949],[114.25663,30.61966]]
df_pos['count'] = 1 ### all stalls that are not multiples have just one sample per stall
df_neg['count'] = 0

## group samples from the same stall
for j in range(len(same_stall)):
    df_pos.loc['comb'+str(j),['x','y']] = same_stall_center[j] #make "combined" stall combining all samples from stall
    if 'A14' in same_stall[j]:
        # df_pos.loc['comb'+str(j),'count'] = 0
        pass
    else:
        df_pos.loc['comb'+str(j),'count'] = len(same_stall[j]) # set count for stalls with multiple positives
    df_pos = df_pos.drop(index=same_stall[j]) #drop the individual stalls in the "combined" samples so we don't double count
df_all = pd.concat((df_pos,df_neg),axis=0)

#live mammals + unk. meat. Concatenated for combined analysis. 
wildlife_stands = gdf[(gdf['group']=='WildlifeVendor') |(gdf['group']=='UnknownMeat')]
unknown_meat = gdf[gdf['group']=='UnknownMeat']

cases= gdf[gdf['group']=='HumanCase']

#get locations of wildlife vendors
series_wildlife = gpd.GeoSeries(wildlife_stands.loc[:,'geometry'])
df_wildlife = pd.DataFrame()
df_wildlife['x'] = series_wildlife.x
df_wildlife['y'] = series_wildlife.y

#get locations of human cases
series_cases = gpd.GeoSeries(cases.loc[:,'geometry'])
df_cases = pd.DataFrame()
df_cases['x'] = series_cases.x
df_cases['y'] = series_cases.y

## calculate distances from each stall to the nearest wildlife vendor/human case
df_all['wild_dist']= [ np.min([distance.distance(tuple(df_all.iloc[j].loc[['y','x']]),tuple(df_wildlife.iloc[k].loc[['y','x']])).meters
                       for k in range(df_wildlife.shape[0])]) for j in range(df_all.shape[0])]
df_all['case_dist']= [ np.min([distance.distance(tuple(df_all.iloc[j].loc[['y','x']]),tuple(df_cases.iloc[k].loc[['y','x']])).meters
                       for k in range(df_cases.shape[0])]) for j in range(df_all.shape[0])]

### Now specify the number of total samples taken per stall

df_pvals = pd.DataFrame()
for N in range(5, 11):
    # Failures is just N - successes
    df_all['fails'] = N-df_all['count']

    ###result in python using statsmodels
    mod1 = smf.glm(formula='count+fails~wild_dist+case_dist', data=df_all, family=sm1.families.Binomial()).fit()
    print('Python Result', mod1.summary())

    ### same test in R with the stats library
    import rpy2.robjects as robjects
    robjects.r("library(stats)")
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri

    # convert to r dataframe
    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_df = robjects.conversion.py2rpy(df_all)

    robjects.globalenv['rdf'] = r_df
    robjects.r('rdf$total = rep('+str(N)+', '+str(df_all.shape[0])+')')

    robjects.globalenv['model'] = robjects.r("glm('cbind(count, total-count)~wild_dist+case_dist', family = binomial, data = rdf)")
    print('binomial model', robjects.r('summary(model)'))
    ## Check deviance to make sure binomial model is justified (i.e. that we don't need quasi-binomial)
    print('p value for nonzero deviance', robjects.r("pchisq(deviance(model),df.residual(model))"))

    vals = list(robjects.r('coef(summary(model))[,4]')[1:])
    coefs = list(np.exp(robjects.r('coef(summary(model))[,1]')[1:]))
    oddsChange = list((1.- np.exp(robjects.r('coef(summary(model))[,1]')[1:]))*100.)
    df_pvals = df_pvals.append(pd.Series(coefs+oddsChange+vals,name='total_'+str(N)))

df_pvals.columns = ['wild_coefs','case_coefs','wild_OddsChange','case_OddsChange','wild_dist_pval','case_dist_pval']
df_pvals.to_csv('pvals_drop2628.csv')

