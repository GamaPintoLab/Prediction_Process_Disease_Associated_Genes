# Thanks to the DoRothEA package, we can easily load the transcription
# factors - target interactions. The output of the load_regulons
# function is an adjacency matrix with the respective HGNC labels
# present in each column and row.


import pandas as pd
import dorothea as dr

pd.set_option('precision', 0)
pd.set_option('display.float_format', lambda x: '%.0f' % x)

dorothea = dr.load_regulons(organism='Human', commercial=False)

# As the interactions come as an adjacency matrix, a conversion to the
# standard output that APID, HuRI and Omnipath present is needed

tf_target_dict = {}
for tf in dorothea.columns:
    tf_targets = list(dorothea[dorothea[tf] == 1][tf].index)
    tf_target_dict[tf] = tf_targets

tf_target_df = pd.DataFrame([(key, var) for (key, L) in tf_target_dict.items() for var in L],
                            columns=['tf', 'target'])

tf_target_df.to_csv('././data/interim/dorothea.csv', header=-1, index=False)
