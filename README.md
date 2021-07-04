Staem5: a novel stacked ensemble method for prediction of m5C site

Before running m5C_predict, users shuold make sure all the following packages are installed in their Python enviroment: 

    numpy == 1.19.5
    pandas == 1.2.4
    sklearn == 0.21.3
    xgboost == 0.90
    lightgbm == 2.3.0
    mlxtend == 0.17.3
    Bio == 0.4.1
    keras==2.30
    tensorflow == 1.14
    python == 3.7

For advanced users who want to perform prediction by using their own data:

 To get the information the user needs to enter for help, run:
 
    python m5C_predict.py --help
    
 or
 
    python m5C_predict.py -h
   
as follows:

>python prom_pred.py -h

Using TensorFlow backend.

usage: m5C_predict.py [-h] --input inputpath [--output OUTPUTFILE]  --species SPECIESFILE  

# Staem5: a novel stacked ensemble method for prediction of m5C site

optional arguments:

  -h, --help                show this help message and exit
  
  --input inputpath     query RNA  sequences to be predicted in fasta format.

  --output OUTPUTFILE       save the prediction results.
  
  --species SPECIESFILE
  
                            --species indicates the specific species, 
                        
                            currently we accept 'Arabidopsis' or 'Mouse'.
                        
  
     


