# Stress Predictor based on Machine Learning


### Learning

With prepared dataset,

```
python learn.py dataset/dataset_learn_v2_005.limma.csv dataset/dataset_learn_v2_005.label.csv --reduce_updown --save output/model_v2 --save_csv output/model_v2 --epoch_count 1000 --batch_size 0
```

* <code>--reduce_updown</code> Use this parameter if dataset is separated into up and down column.
* <code>--save</code> Directory to save model data
* <code>--save_csv</code> Directory to save model data (in csv form)
* <code>--epoch_count</code> epoch to learn
* <code>--batch_size</code> learning batch size

### Testing prediction /w learnt parameter

```
python test.py dataset/dataset_learn_v2_005.limma.csv dataset/dataset_learn_v2_005.label.csv --load output/model_v2 --save output/model_v2_test
python test.py dataset/dataset_test_005.limma.csv dataset/dataset_test_005.label.csv --load output/model_v2 --save output/model_v2_test
```

* <code>--load</code> Directory of learn.py output
* <code>--save</code> Directory to save test result

### GSEA (Gene set enrichment test) Testing

```
python GOenrichment.py $1/output/params_generator.csv --trait2genes GOBPname2gene.arabidopsis.txt --column_name "heat,salt,drought,cold" --count_cut 500 -o $1/output/GSEA_learn_top500.csv --descending --max_trait_cut 1
```

* <code>--trait2genes</code> Path to trait-to-gene file.
* <code>--column_name</code> Column names including p-value.
* <code>--count_cut</code> How many gene set to be shown?
* <code>-o</code> Path to output
* <code>--descending</code> Result data in descending p-value order
* <code>--max-trait-cut</code> Remove genes with too little traits.


### How to generate dataset?

Just prepare your expression matrix (CEL file processed) and metadata (may write TSD header or generate from CSV). then use these commands:

```
# split total dataset file into time-series-data formatted file
python ../dataset.py tsd -f cold.csv
python ../dataset.py tsd -f heat.csv
# read directory just for test
python ../dataset.py read -d ./

# without any processing, just output raw p-value
python ../dataset.py compile -d ./ --tool limma -t dataset_test_005
# processing p-value /w threshold .05, (signed)
python ../dataset.py compile -d ./ --tool limma -t dataset_test_005 --pvalue 0.05
# processing p-value /w threshold .05, with reindexing & separating up/down signal (unsigned)
python ../dataset.py compile -d ./ --tool limma -t dataset_test_005 --pvalue 0.05 --updown --reindex_df ../dataprocess/ttest_pval.csv

# generate gene list order for comparison with other result (optional)
python dataset.py gen_genelist -f data/GSE3326_1.tsd -t dataset/index
# create dataset with pvalue 0.05
python dataset.py compile -d data_test/ --tool limma -t dataset/dataset_test_005 --pvalue 0.05 --reindex_file dataset/index.txt --updown --label_index "heat,salt,drought,cold"
# create dataset with additional filter
python dataset.py compile -d data/ --tool limma -t dataset/dataset_learn_v2_005 --pvalue 0.05 --reindex_file dataset/index.txt  --updown --label_index "heat,salt,drought,cold" --filters "Species:Arabidopsis Thaliana,MinRepCnt:2"
```

Then DEG infomation will be stored at `dataset_test.csv` matrix.
You can use `--filter` to select TSData with specific condition (with column name and value).
