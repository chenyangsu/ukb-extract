# ukb-extract

## Usage
        1. Open extract.py and set paths and user-options
        2. python extract.py

Data-field names in UKB file are coded as
<pre>f . field-id . visit number . instance number</pre>
<pre>ex) f.34.0.0</pre>

They can be one of the following types:
* Integer
* Categorical (single)
* Categorical (multiple)
* Continuous
* Text
* Date
* Time
* Compound

## <b>extract.py</b>
1. Webscrape descriptions, categories and variable-types for each field-id from: 
2. http://biobank.ctsu.ox.ac.uk/crystal/list.cgi.
3. Enumerate field-id’s for each variable type.
    <pre>
    * Output/Integer.datafields
    * Output/Categorical_single.datafields
    * Output/Categorical_multiple.datafields
    * Output/Continuous.datafields
    </pre>
4. Split UK Biobank by variable types: Integer, Categorical and Continuous. 
    * For continuous and integer features, encode missing values as ‘NA’.
        <pre>
        * Output/Integer.phenoslice.csv
        * Output/Continuous.phenoslice.csv
        </pre>
    * For categorical features, 
        * Enumerate all possible categories.
        * One-hot encode all field-id/category combinations. 
        <pre>
            * Output/Categorical_single.phenoslice.csv
            * Output/Categorical_multiple.phenoslice.csv
        </pre>
5. Merge dataset into a single file.
    * Output/features_merged.csv

### Note:
* ‘f.xxxxx.x.x’ represents integer or continuous columns
    * ex) f.3064.0.1 would represent the measurement #1 of data-field 3064 on assessment day 0.
* ‘c.xxxx.x’ represents categorical columns.
    * ex) c.3606.1 would be the indicator that an individual is in category 1 (‘Yes’) for data-field 3606.

To look up the significance of a field, use the UK Biobank data showcase: http://biobank.ndph.ox.ac.uk/showcase/
