from flask import Flask, render_template, request, flash, redirect, url_for
from joblib import load
from utils.EDA_utils import *
from utils.utils import *
import os
from werkzeug.utils import secure_filename
import vcf
import subprocess
UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = {'vcf', 'gz'}

app= Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = b'_5#y2L"F4Q8z\n\xec]/'
model=load('model.joblib')

def clean_cols(df):
    model_cols=set([x[:-7] if '_is_nan' in x else x for x in model.feature_names_in_])
    missing_cols=model_cols-set(df.columns)
    for col in missing_cols:
        df[col]=pd.Series(dtype=float)
    
    df=df.loc[:, model_cols]

    return df
    
def fix_nan(df):
    missing_nan=set(model.feature_names_in_)-set(df.columns)
    for x in missing_nan:
        print(x)
    for nan in missing_nan:
        df[nan]=0
    return df.loc[:, model.feature_names_in_]
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def clean_files(folder):
    for root, dirs, files in os.walk(folder):
        for file in files:
            # Remove sample files
            if "sample" in file:
                os.remove(folder / file)
            # Rewrite genotype info without csq
            elif "genotype" in file:
                df = pd.read_csv(folder / file, sep="\t")
                df.drop("CSQ", axis=1).to_csv(folder / f"{file[:6]}_genotype.csv.gz", sep=";", compression="gzip")
                os.remove(folder / file)
            # Change default separator and compress
            elif "default" in file:
                df = pd.read_csv(folder / file, sep="\t")
                df.to_csv(folder / f"{file[:6]}_default.csv.gz", sep=";", compression="gzip")
                os.remove(folder / file)

def move_csq(file):
    data = []

    vcf_file = file
    with open(vcf_file, 'rb') as vcf_file_binary:
        vcf_reader = vcf.Reader(vcf_file_binary)
        
        column_csq_headers = vcf_reader.infos["CSQ"].desc[50:].split("|")

        for record in vcf_reader:
            for field in vcf_reader.infos.keys():
                if field == "CSQ":
                    data.append(record.INFO.get(field)[0].split("|"))
                    break

    df = pd.DataFrame(data, columns=column_csq_headers)
    df.to_csv(f"{file.split('.')[0]}_csq.csv.gz", sep=";", compression="gzip")


@app.route('/')
def index():
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload():
    if 'file' not in request.files:
        flash('No file part')
        return redirect(url_for('index'))

    file = request.files['file']

    if file.filename == '':
        flash('No selected file')
        return redirect(url_for('index'))

    
    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        print(filename.split('.')[0])
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(file_path)
        command=f'unvcf {file_path} {app.config["UPLOAD_FOLDER"]}'
        subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        clean_files(Path(app.config["UPLOAD_FOLDER"]))
        move_csq(file_path)
        data=get_dataset(app.config['UPLOAD_FOLDER'], [filename.split('.')[0]], 'both', option_csq="potential", 
            options_genotype=["potential", "all"], with_default=True)
        data=clean_cols(data)
        data=preprocess(data)
        data=fix_nan(data)
        data, encoders1, target_mapping1 = encode(data, 'ACMG_class', {'Benign': 0,'Likely%40Benign': 1,'Uncertain%40Significance': 2,'Likely%40Pathogenic': 3,'Pathogenic': 4})
        for x, y in zip(sorted(model.feature_names_in_), sorted(data.columns)):
            print(x,y, x==y, data[y].isnull().any())
        y_pred=model.predict(data)
        return render_template('results.html', results=y_pred)
    else:
        flash('Invalid file type')
        return redirect(url_for('index'))

    
if __name__ == '__main__':
    app.run(debug=True)