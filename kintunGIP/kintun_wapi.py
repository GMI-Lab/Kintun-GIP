# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 16:24:44 2021

@author: roberto
"""

import sys
import subprocess
import csv
from Bio import SeqIO
from Bio import SeqFeature as sf
from Bio.SeqRecord import SeqRecord
import glob
import os
import datetime
import json

# Primera función, chequear las coordenadas



# Nombre del archivo multifasta de entrada del programa
in_file = sys.argv[1]
position_in_file_name = in_file.rfind(".")
in_file_short = in_file[:position_in_file_name]

# Crear carpeta donde estarán todos los resultados
path_file = in_file_short
try:
    os.mkdir(path_file)
except OSError:
    print("Directory /%s not created, because already exists." % path_file)
else:
    print("Successfully created the directory /%s." % path_file)
print(' ')

##############################################################################
############# 1. GENERACIÓN DE ARCHIVOS SEPARADOS PARA CADA tDNA #############
##############################################################################

# En cada paso del programa va a ir una línea con el nombre del paso y la hora exacta

print(' ')
actual_time = datetime.datetime.now().strftime('%H:%M:%S')
print('[' + actual_time + '] ===== Paso 1: Creación de archivos para cada tDNA y no redundantes. INICIANDO ======')
print(' ')

# 1.1 Crear una lista con los tDNA presentes en el multifasta

all_sites = []
with open(in_file, 'r') as opened_file:
    for record in SeqIO.parse(opened_file, "fasta"):
        ident = record.id
        all_sites.append(ident.split("-")[0])

# 1.2 Hacer una lista no redundante de nombres de tDNA

non_redundant_sites = []
for site in all_sites:
    if site not in non_redundant_sites:
        non_redundant_sites.append(site)
    else:
        continue
'''
# 1.3 Crear la carpeta que contendrá las islas agrupadas por GI

path_separ = "tDNAs_GIs_separately"
path = os.path.join(path_file, path_separ)
try:
    os.mkdir(path)
except OSError:
    print ("Directory %s not created, because already exists." % path)
else:
    print ("Successfully created the directory %s." % path)

# 1.4 Crear multifastas diferentes para cada tDNA dentro de la carpeta recién creada

islands = []
for site in non_redundant_sites:
    site_islands = []
    with open(in_file, 'r') as content:
        for record in SeqIO.parse(content, "fasta"):
            if site in record.id:
                site_islands.append(record)
                islands.append(record.id)
            else:
                continue
    for island in site_islands:
        with open("./%s/tDNAs_GIs_separately/%s_GIs.fasta" %(in_file_short, site), 'a') as output_content:
            SeqIO.write(island,output_content,"fasta")


# 1.5 Crear o sobreescribir archivo .txt que tendrá un resumen del procesamiento del archivo

with open('./%s/%s_final_annotation.txt' %(in_file_short,in_file_short), 'w') as txt_output:
    txt_output.write("")

# 1.5.1 Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' %(in_file_short,in_file_short), 'a') as txt_output:
    txt_output.write("Summary of the analysis of: " + in_file + "\n\n")

# 1.5.2 Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' %(in_file_short,in_file_short), 'a') as txt_output:
    txt_output.write("1. The file has " + str(len(islands)) + " genomic islands inserted in " + str(len(non_redundant_sites)) + " different tDNAs.")


# Fin de la primera parte

print(' ')
actual_time = datetime.datetime.now().strftime('%H:%M:%S')
print('[' + actual_time + '] ===== Paso 1: Creación de archivos para cada tDNA y no redundantes. FINALIZADO =====')
print(' ')



##############################################################################
################# 2. EJECUTAR PROKKA CON LAS BASES DE DATOS ##################
##############################################################################

# Inicio de la segunda parte

print(' ')
actual_time = datetime.datetime.now().strftime('%H:%M:%S')
print('[' + actual_time + '] ======== Paso 2: Anotación general con diferentes bases de datos. INICIANDO ========')
print(' ')

# 2.1 EJECUCIÓN: Tomar cada multifasta, ejecutar PROKKA junto a las bases de datos

for site in non_redundant_sites:
    cmd_prokka_1 = "prokka --proteins db_VFDB_prokka.fasta --outdir ./%s/results_prokka/%s_GIs/VFDB --prefix %s_VFDB %s/%s_GIs.fasta" %(in_file_short,site,site,path,site)
    ex_prokka_1 = subprocess.run(cmd_prokka_1, shell=True)

    cmd_prokka_2 = "prokka --proteins db_integrase_ncbi_prokka.fasta --outdir ./%s/results_prokka/%s_GIs/integrase --prefix %s_integrase %s/%s_GIs.fasta" %(in_file_short,site,site,path,site)
    ex_prokka_2 = subprocess.run(cmd_prokka_2, shell=True)

    cmd_prokka_3 = "prokka --proteins db_CARD_prokka.fasta --outdir ./%s/results_prokka/%s_GIs/CARD --prefix %s_CARD %s/%s_GIs.fasta" %(in_file_short,site,site,path,site)
    ex_prokka_3 = subprocess.run(cmd_prokka_3, shell=True)

    cmd_prokka_4 = "prokka --proteins db_TADB_prokka.fasta --outdir ./%s/results_prokka/%s_GIs/TADB --prefix %s_TADB %s/%s_GIs.fasta" %(in_file_short,site,site,path,site)
    ex_prokka_4 = subprocess.run(cmd_prokka_4, shell=True)

    cmd_prokka_5 = "prokka --proteins db_BacMet2_prokka.fasta --outdir ./%s/results_prokka/%s_GIs/BacMet2 --prefix %s_BacMet2 %s/%s_GIs.fasta" %(in_file_short,site,site,path,site)
    ex_prokka_5 = subprocess.run(cmd_prokka_5, shell=True)
'''

# Aprovechando el análisis reciente, contar la cantidad de CDS y tRNA
c_num_CDS = 0
c_num_tRNA = 0

# Para contar cuántos CDS se analizarán en las siguientes etapas
for site in non_redundant_sites:
    # Archivo del que se obtendrá la información, da lo mismo de cuál base de datos
    in_CDS_gbk = "./%s/results_prokka/%s_GIs/VFDB/%s_VFDB.gbk" % (in_file_short, site, site)
    # Analizar los archivos de cada carpeta
    with open(in_CDS_gbk, "r") as annot_file:
        # Se le pide a SeqIO que manipule la info del archivo como un genbank
        for record in SeqIO.parse(annot_file, "genbank"):
            for rf in record.features:
                if rf.type == "CDS":
                    c_num_CDS = c_num_CDS + 1
                elif rf.type == "tRNA":
                    c_num_tRNA = c_num_tRNA + 1
                else:
                    continue

# Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' % (in_file_short, in_file_short), 'a') as txt_output:
    txt_output.write(
        " Inside the genomic islands, " + str(c_num_CDS) + " CDS and " + str(c_num_tRNA) + " tRNAs were found.\n\n")

# 2.2 PROCESAMIENTO: Dejar el output de PROKKA en un formato que pueda ser
#     tomado por otro comando y lo agregue al resultado

# Creación de carpeta que contendrá todos los outputs de los procesamientos

path_outputs = "outputs"
path = os.path.join(path_file, path_outputs)
try:
    os.mkdir(path)
except OSError:
    print("Directory /%s not created, because already exists." % path)
else:
    print("Successfully created the directory /%s." % path)
print(' ')

# 2.2.1 Creación del archivo output para datos de integrase NCBI

# Crear el archivo nuevo, sobreescribiendo si ya existía
header_prokka_tsv = ["tDNA", "Start", "End", "Strand", "ID", "Inference", "Product"]
with open("./%s/outputs/output_integrase.tsv" % in_file_short, 'w', newline='') as f_output:
    tsv_output = csv.writer(f_output, delimiter='\t')
    tsv_output.writerow(header_prokka_tsv)

# Analizar los archivos de cada carpeta
for site in non_redundant_sites:
    prokka_in_file = "./%s/results_prokka/%s_GIs/integrase/%s_integrase.gff" % (in_file_short, site, site)
    # Lista que contendrá toda la información de integrase NCBI y que será escrita en el output
    list_of_genes = []

    # Análisis de cada archivo .gff
    with open(prokka_in_file, 'r') as data_gff:
        for line in data_gff:

            if '##FASTA' in line:
                # Está al final del archivo, por lo tanto, se termina el análisis
                break
            elif '##' in line:
                # Se salta esta línea
                continue

            else:
                # Se procesa la información de esta línea
                gene_prokka = line.split("\t")
                last_column = gene_prokka[8]
                last_column_sp = last_column.split(";")
                attributes = {}
                for item in last_column_sp:
                    att = item.split("=")

                    # Eliminar el salto de línea
                    if att[0] == "product":
                        list_pre = att[1]
                        list_post = list_pre[:-1]
                        att[1] = list_post
                    attributes[att[0]] = att[1]

                # Se junta toda la información útil en una variable
                useful_info = [gene_prokka[0], gene_prokka[3], gene_prokka[4], gene_prokka[6], attributes["ID"],
                               attributes["inference"], attributes["product"]]
                # Solo los identificados por esta base de datos, serán almacenados
                if "db_integrase_ncbi_prokka.fasta" in attributes["inference"]:
                    list_of_genes.append(useful_info)
                else:
                    continue

    # Escribir la información en el archivo output
    with open("./%s/outputs/output_integrase.tsv" % in_file_short, 'a+', newline='') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        tsv_output.writerows(list_of_genes)

    # 2.2.2 Creación del archivo output para datos de TADB

# Crear el archivo nuevo, sobreescribiendo si ya existía
header_prokka_tsv = ["tDNA", "Start", "End", "Strand", "ID", "Inference", "Product"]
with open("./%s/outputs/output_TADB.tsv" % in_file_short, 'w', newline='') as f_output:
    tsv_output = csv.writer(f_output, delimiter='\t')
    tsv_output.writerow(header_prokka_tsv)

# Analizar los archivos de cada carpeta
for site in non_redundant_sites:
    prokka_in_file = "./%s/results_prokka/%s_GIs/TADB/%s_TADB.gff" % (in_file_short, site, site)
    # Lista que contendrá toda la información de TADB y que será escrita en el output
    list_of_genes = []

    # Análisis de cada archivo .gff
    with open(prokka_in_file, 'r') as data_gff:
        for line in data_gff:

            if '##FASTA' in line:
                # Está al final del archivo, por lo tanto, se termina el análisis
                break
            elif '##' in line:
                # Se salta esta línea
                continue

            else:
                # Se procesa la información de esta línea
                gene_prokka = line.split("\t")
                last_column = gene_prokka[8]
                last_column_sp = last_column.split(";")
                attributes = {}
                for item in last_column_sp:
                    att = item.split("=")

                    # Eliminar el salto de línea
                    if att[0] == "product":
                        list_pre = att[1]
                        list_post = list_pre[:-1]
                        att[1] = list_post
                    attributes[att[0]] = att[1]

                # Se junta toda la información útil en una variable
                useful_info = [gene_prokka[0], gene_prokka[3], gene_prokka[4], gene_prokka[6], attributes["ID"],
                               attributes["inference"], attributes["product"]]
                # Solo los identificados por esta base de datos, serán almacenados
                if "db_TADB_prokka.fasta" in attributes["inference"]:
                    list_of_genes.append(useful_info)
                else:
                    continue

    # Escribir la información en el archivo output
    with open("./%s/outputs/output_TADB.tsv" % in_file_short, 'a+', newline='') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        tsv_output.writerows(list_of_genes)

    # 2.2.3 Creación del archivo output para datos de VFDB

# Crear el archivo nuevo, sobreescribiendo si ya existía
header_prokka_tsv = ["tDNA", "Start", "End", "Strand", "ID", "Inference", "Product"]
with open("./%s/outputs/output_VFDB.tsv" % in_file_short, 'w', newline='') as f_output:
    tsv_output = csv.writer(f_output, delimiter='\t')
    tsv_output.writerow(header_prokka_tsv)

# Analizar los archivos de cada carpeta
for site in non_redundant_sites:
    prokka_in_file = "./%s/results_prokka/%s_GIs/VFDB/%s_VFDB.gff" % (in_file_short, site, site)
    # Lista que contendrá toda la información de VFDB y que será escrita en el output
    list_of_genes = []

    # Análisis de cada archivo .gff
    with open(prokka_in_file, 'r') as data_gff:
        for line in data_gff:

            if '##FASTA' in line:
                # Está al final del archivo, por lo tanto, se termina el análisis
                break
            elif '##' in line:
                # Se salta esta línea
                continue

            else:
                # Se procesa la información de esta línea
                gene_prokka = line.split("\t")
                last_column = gene_prokka[8]
                last_column_sp = last_column.split(";")
                attributes = {}
                for item in last_column_sp:
                    att = item.split("=")

                    # Eliminar el salto de línea
                    if att[0] == "product":
                        list_pre = att[1]
                        list_post = list_pre[:-1]
                        att[1] = list_post
                    attributes[att[0]] = att[1]

                # Se junta toda la información útil en una variable
                useful_info = [gene_prokka[0], gene_prokka[3], gene_prokka[4], gene_prokka[6], attributes["ID"],
                               attributes["inference"], attributes["product"]]
                # Solo los identificados por esta base de datos, serán almacenados
                if "db_VFDB_prokka.fasta" in attributes["inference"]:
                    list_of_genes.append(useful_info)
                else:
                    continue

    # Escribir la información en el archivo output
    with open("./%s/outputs/output_VFDB.tsv" % in_file_short, 'a+', newline='') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        tsv_output.writerows(list_of_genes)

    # 2.2.4 Creación del archivo output para datos de CARD

# Crear el archivo nuevo, sobreescribiendo si ya existía
header_prokka_tsv = ["tDNA", "Start", "End", "Strand", "ID", "Inference", "Product"]
with open("./%s/outputs/output_CARD.tsv" % in_file_short, 'w', newline='') as f_output:
    tsv_output = csv.writer(f_output, delimiter='\t')
    tsv_output.writerow(header_prokka_tsv)

# Analizar los archivos de cada carpeta
for site in non_redundant_sites:
    prokka_in_file = "./%s/results_prokka/%s_GIs/CARD/%s_CARD.gff" % (in_file_short, site, site)
    # Lista que contendrá toda la información de CARD y que será escrita en el output
    list_of_genes = []

    # Análisis de cada archivo .gff
    with open(prokka_in_file, 'r') as data_gff:
        for line in data_gff:

            if '##FASTA' in line:
                # Está al final del archivo, por lo tanto, se termina el análisis
                break
            elif '##' in line:
                # Se salta esta línea
                continue

            else:
                # Se procesa la información de esta línea
                gene_prokka = line.split("\t")
                last_column = gene_prokka[8]
                last_column_sp = last_column.split(";")
                attributes = {}
                for item in last_column_sp:
                    att = item.split("=")

                    # Eliminar el salto de línea
                    if att[0] == "product":
                        list_pre = att[1]
                        list_post = list_pre[:-1]
                        att[1] = list_post
                    attributes[att[0]] = att[1]

                # Se junta toda la información útil en una variable
                useful_info = [gene_prokka[0], gene_prokka[3], gene_prokka[4], gene_prokka[6], attributes["ID"],
                               attributes["inference"], attributes["product"]]
                # Solo los identificados por esta base de datos, serán almacenados
                if "db_CARD_prokka.fasta" in attributes["inference"]:
                    list_of_genes.append(useful_info)
                else:
                    continue

                    # Escribir la información en el archivo output
    with open("./%s/outputs/output_CARD.tsv" % in_file_short, 'a+', newline='') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        tsv_output.writerows(list_of_genes)

    # 2.2.5 Creación del archivo output para datos de BacMet2

# Crear el archivo nuevo, sobreescribiendo si ya existía
header_prokka_tsv = ["tDNA", "Start", "End", "Strand", "ID", "Inference", "Product"]
with open("./%s/outputs/output_BacMet2.tsv" % in_file_short, 'w', newline='') as f_output:
    tsv_output = csv.writer(f_output, delimiter='\t')
    tsv_output.writerow(header_prokka_tsv)

# Analizar los archivos de cada carpeta
for site in non_redundant_sites:
    prokka_in_file = "./%s/results_prokka/%s_GIs/BacMet2/%s_BacMet2.gff" % (in_file_short, site, site)
    # Lista que contendrá toda la información de BacMet2 y que será escrita en el output
    list_of_genes = []

    # Análisis de cada archivo .gff
    with open(prokka_in_file, 'r') as data_gff:
        for line in data_gff:

            if '##FASTA' in line:
                # Está al final del archivo, por lo tanto, se termina el análisis
                break
            elif '##' in line:
                # Se salta esta línea
                continue

            else:
                # Se procesa la información de esta línea
                gene_prokka = line.split("\t")
                last_column = gene_prokka[8]
                last_column_sp = last_column.split(";")
                attributes = {}
                for item in last_column_sp:
                    att = item.split("=")

                    # Eliminar el salto de línea
                    if att[0] == "product":
                        list_pre = att[1]
                        list_post = list_pre[:-1]
                        att[1] = list_post
                    attributes[att[0]] = att[1]

                # Se junta toda la información útil en una variable
                useful_info = [gene_prokka[0], gene_prokka[3], gene_prokka[4], gene_prokka[6], attributes["ID"],
                               attributes["inference"], attributes["product"]]
                # Solo los identificados por esta base de datos, serán almacenados
                if "db_BacMet2_prokka.fasta" in attributes["inference"]:
                    list_of_genes.append(useful_info)
                else:
                    continue

    # Escribir la información en el archivo output
    with open("./%s/outputs/output_BacMet2.tsv" % in_file_short, 'a+', newline='') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        tsv_output.writerows(list_of_genes)

# Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' % (in_file_short, in_file_short), 'a') as txt_output:
    txt_output.write("2. ANNOTATION: PROKKA analysis used several databases:\n")

# 2.3 INCORPORACIÓN: agregar toda la información de los outputs en un archivo
#     genbank para que recopile todas las anotaciones de PROKKA, según la
#     importancia de cada base de datos


# 2.3.1 Se incorpora en el qualifier /note los CDS anotados por PROKKA para
#       la primera base de datos que es integrase NCBI y aunque cambie la anotación
#       más adelante, quede registro que fue también anotado con integrase NCBI.

# Contador de CDS incluidos a los datos (c de count)
c_cds = 0

# Se crea un loop para que analice los archivos .gbk que servirán como base,
# para después modificarlos con la información del output correspondiente.
for site in non_redundant_sites:
    # Archivo del cual se obtendrá la información usada como base para la anotación
    in_prokka_gbk = "./%s/results_prokka/%s_GIs/integrase/%s_integrase.gbk" % (in_file_short, site, site)
    # Archivo con las anotaciones tabuladas que fueron realizadas con integrase NCBI
    in_data = "./%s/outputs/output_integrase.tsv" % (in_file_short)

    # Crear el archivo nuevo, sobreescribiendo si ya existía
    record = ""
    out_prokka_gbk = "./%s/results_prokka/%s_GIs/replace_0_from_integrase_to_integrase.gb" % (in_file_short, site)
    with open(out_prokka_gbk, "w") as gbk_output:
        SeqIO.write(record, gbk_output, "genbank")

    # Analizar los archivos de cada carpeta
    with open(in_prokka_gbk, "r") as annot_file:
        # Se le pide a SeqIO que manipule la info del archivo como un genbank
        for record in SeqIO.parse(annot_file, "genbank"):
            # El objeto record es la secuencia y las anotaciones que ya tiene
            # La sequencia está en record.seq y las anotaciones en record.features

            # Abrir la tabla de las anotaciones
            with open(in_data, "r") as data_genes:
                for line in data_genes:
                    sep = line.strip().split("\t")

                    if sep[0] == "tDNA":
                        # Se salta esta línea
                        continue

                    else:
                        # Se procesa la información de esta línea
                        gi = sep[0]
                        start = int(sep[1])
                        end = int(sep[2])
                        # Se coloca la hebra en el formato adecuado
                        if sep[3] == "+":
                            sense = +1
                        else:
                            sense = -1

                    # Se encuentra la correspondencia entre la isla del output y del .gbk
                    if gi == record.id:
                        # Se crea una lista con las coordenadas del CDS
                        spos = sf.ExactPosition(start - 1)
                        epos = sf.ExactPosition(end)
                        l1 = sf.FeatureLocation(spos, epos, strand=sense)

                        # Se lee cada entrada del genbank
                        for rf in record.features:
                            # Asegurarse que el CDS coincida con el del output
                            if l1 == rf.location:
                                # Probar si existe el qualifier /note, porque
                                # el resto de la anotación ya está incorporada
                                try:
                                    # Modificar cada diccionario para agregar solamente el /note
                                    list_note = rf.qualifiers['note']
                                    list_note.append("identified in NCBI integrases database")
                                    rf.qualifiers['note'] = list_note
                                    c_cds = c_cds + 1
                                except KeyError:
                                    # Si no existe /note, se crea
                                    rf.qualifiers['note'] = "identified in NCBI integrases database"
                                    c_cds = c_cds + 1
                            else:
                                continue

            # Escribir la información en el archivo gbk que agrupa toda la información
            with open(out_prokka_gbk, "a") as gbk_output:
                SeqIO.write(record, gbk_output, "genbank")

# Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' % (in_file_short, in_file_short), 'a') as txt_output:
    txt_output.write("   " + str(c_cds) + " CDS annotated with NCBI integrases database.\n")

print(str(c_cds) + " CDS annotated with NCBI integrases database were included in the data.")

# 2.3.2 Se incorpora en el qualifier /note los CDS anotados por PROKKA para
#       la segunda base de datos que es TADB y aunque cambie la anotación
#       más adelante, quede registro que fue también anotado con TADB.

# Contador de CDS incluidos a los datos
c_cds = 0

# En este caso, de la primera modificación usaremos como base los archivos integrase NCBI
# y le sumaremos la anotación de los resultados TADB.
for site in non_redundant_sites:
    # Archivo del cual se obtendrá la información usada como base para la anotación
    in_prokka_gbk = "./%s/results_prokka/%s_GIs/replace_0_from_integrase_to_integrase.gb" % (in_file_short, site)
    # Archivo con las anotaciones tabuladas que fueron realizadas con TADB
    in_data = "./%s/outputs/output_TADB.tsv" % (in_file_short)

    # Crear el archivo nuevo, sobreescribiendo si ya existía
    record = ""
    out_prokka_gbk = "./%s/results_prokka/%s_GIs/replace_1_from_integrase_to_TADB.gb" % (in_file_short, site)
    with open(out_prokka_gbk, "w") as gbk_output:
        SeqIO.write(record, gbk_output, "genbank")

    # Analizar los archivos de cada carpeta
    with open(in_prokka_gbk, "r") as annot_file:
        # Se le pide a SeqIO que manipule la info del archivo como un genbank
        for record in SeqIO.parse(annot_file, "genbank"):
            # El objeto record es la secuencia y las anotaciones que ya tiene
            # La sequencia está en record.seq y las anotaciones en record.features

            # Abrir la tabla de las anotaciones
            with open(in_data, "r") as data_genes:
                for line in data_genes:
                    sep = line.strip().split("\t")

                    if sep[0] == "tDNA":
                        # Se salta esta línea
                        continue

                    else:
                        # Se procesa la información de esta línea
                        gi = sep[0]
                        start = int(sep[1])
                        end = int(sep[2])
                        # Se coloca la hebra en el formato adecuado
                        if sep[3] == "+":
                            sense = +1
                        else:
                            sense = -1

                    # Se encuentra la correspondencia entre la isla del output y del .gbk
                    if gi == record.id:
                        # Se crea una lista con las coordenadas del CDS
                        spos = sf.ExactPosition(start - 1)
                        epos = sf.ExactPosition(end)
                        l1 = sf.FeatureLocation(spos, epos, strand=sense)

                        # Se lee cada entrada del genbank
                        for rf in record.features:

                            # Asegurarse que el CDS coincida con el del output
                            if l1 == rf.location:
                                # Modificar cada diccionario para agregar inference, product y note
                                pre_inference_list = sep[5]
                                inference_list = list(pre_inference_list.split(','))
                                rf.qualifiers['inference'] = inference_list
                                rf.qualifiers['product'] = sep[6]

                                # Probar si existe el qualifier /note,
                                try:
                                    # Modificar cada diccionario para agregar solamente el /note
                                    list_note = rf.qualifiers['note']
                                    list_note.append("identified in TADB database")
                                    rf.qualifiers['note'] = list_note
                                    c_cds = c_cds + 1
                                except KeyError:
                                    # Si no existe /note, se crea
                                    rf.qualifiers['note'] = "identified in TADB database"
                                    c_cds = c_cds + 1
                            else:
                                continue

            # Escribir la información en el archivo gbk que agrupa toda la información
            with open(out_prokka_gbk, "a") as gbk_output:
                SeqIO.write(record, gbk_output, "genbank")

# Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' % (in_file_short, in_file_short), 'a') as txt_output:
    txt_output.write("   " + str(c_cds) + " CDS annotated with TADB database.\n")

print(str(c_cds) + " CDS annotated with TADB database were included in the data.")

# 2.3.3 Se incorpora en el qualifier /note los CDS anotados por PROKKA para
#       la tercera base de datos que es VFDB y aunque cambie la anotación
#       más adelante, quede registro que fue también anotado con VFDB.

# Contador de CDS incluidos a los datos
c_cds = 0

# En este caso, de la segunda modificación usaremos como base los archivos TADB
# y le sumaremos la anotación de los resultados VFDB.
for site in non_redundant_sites:
    # Archivo del cual se obtendrá la información usada como base para la anotación
    in_prokka_gbk = "./%s/results_prokka/%s_GIs/replace_1_from_integrase_to_TADB.gb" % (in_file_short, site)
    # Archivo con las anotaciones tabuladas que fueron realizadas con VFDB
    in_data = "./%s/outputs/output_VFDB.tsv" % (in_file_short)

    # Crear el archivo nuevo, sobreescribiendo si ya existía
    record = ""
    out_prokka_gbk = "./%s/results_prokka/%s_GIs/replace_2_from_TADB_to_VFDB.gb" % (in_file_short, site)
    with open(out_prokka_gbk, "w") as gbk_output:
        SeqIO.write(record, gbk_output, "genbank")

    # Analizar los archivos de cada carpeta
    with open(in_prokka_gbk, "r") as annot_file:
        # Se le pide a SeqIO que manipule la info del archivo como un genbank
        for record in SeqIO.parse(annot_file, "genbank"):
            # El objeto record es la secuencia y las anotaciones que ya tiene
            # La sequencia está en record.seq y las anotaciones en record.features

            # Abrir la tabla de las anotaciones
            with open(in_data, "r") as data_genes:
                for line in data_genes:
                    sep = line.strip().split("\t")

                    if sep[0] == "tDNA":
                        # Se salta esta línea
                        continue

                    else:
                        # Se procesa la información de esta línea
                        gi = sep[0]
                        start = int(sep[1])
                        end = int(sep[2])
                        # Se coloca la hebra en el formato adecuado
                        if sep[3] == "+":
                            sense = +1
                        else:
                            sense = -1

                    # Se encuentra la correspondencia entre la isla del output y del .gbk
                    if gi == record.id:
                        # Se crea una lista con las coordenadas del CDS
                        spos = sf.ExactPosition(start - 1)
                        epos = sf.ExactPosition(end)
                        l1 = sf.FeatureLocation(spos, epos, strand=sense)

                        # Se lee cada entrada del genbank
                        for rf in record.features:

                            # Asegurarse que el CDS coincida con el del output
                            if l1 == rf.location:
                                # Modificar cada diccionario para agregar inference, product y note
                                pre_inference_list = sep[5]
                                inference_list = list(pre_inference_list.split(','))
                                rf.qualifiers['inference'] = inference_list
                                rf.qualifiers['product'] = sep[6]

                                # Probar si existe el qualifier /note,
                                try:
                                    # Modificar cada diccionario para agregar solamente el /note
                                    list_note = rf.qualifiers['note']
                                    list_note.append("identified in VFDB database")
                                    rf.qualifiers['note'] = list_note
                                    c_cds = c_cds + 1
                                except KeyError:
                                    # Si no existe /note, se crea
                                    rf.qualifiers['note'] = "identified in VFDB database"
                                    c_cds = c_cds + 1
                            else:
                                continue

            # Escribir la información en el archivo gbk que agrupa toda la información
            with open(out_prokka_gbk, "a") as gbk_output:
                SeqIO.write(record, gbk_output, "genbank")

# Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' % (in_file_short, in_file_short), 'a') as txt_output:
    txt_output.write("   " + str(c_cds) + " CDS annotated with VFDB database.\n")

print(str(c_cds) + " CDS annotated with VFDB database were included in the data.")

# 2.3.4 Se incorpora en el qualifier /note los CDS anotados por PROKKA para
#       la cuarta base de datos que es CARD y aunque cambie la anotación
#       más adelante, quede registro que fue también anotado con CARD.

# Contador de CDS incluidos a los datos
c_cds = 0

# En este caso, de la tercera modificación usaremos como base los archivos VFDB
# y le sumaremos la anotación de los resultados CARD.
for site in non_redundant_sites:
    # Archivo del cual se obtendrá la información usada como base para la anotación
    in_prokka_gbk = "./%s/results_prokka/%s_GIs/replace_2_from_TADB_to_VFDB.gb" % (in_file_short, site)
    # Archivo con las anotaciones tabuladas que fueron realizadas con CARD
    in_data = "./%s/outputs/output_CARD.tsv" % (in_file_short)

    # Crear el archivo nuevo, sobreescribiendo si ya existía
    record = ""
    out_prokka_gbk = "./%s/results_prokka/%s_GIs/replace_3_from_VFDB_to_CARD.gb" % (in_file_short, site)
    with open(out_prokka_gbk, "w") as gbk_output:
        SeqIO.write(record, gbk_output, "genbank")

    # Analizar los archivos de cada carpeta
    with open(in_prokka_gbk, "r") as annot_file:
        # Se le pide a SeqIO que manipule la info del archivo como un genbank
        for record in SeqIO.parse(annot_file, "genbank"):
            # El objeto record es la secuencia y las anotaciones que ya tiene
            # La sequencia está en record.seq y las anotaciones en record.features

            # Abrir la tabla de las anotaciones
            with open(in_data, "r") as data_genes:
                for line in data_genes:
                    sep = line.strip().split("\t")

                    if sep[0] == "tDNA":
                        # Se salta esta línea
                        continue

                    else:
                        # Se procesa la información de esta línea
                        gi = sep[0]
                        start = int(sep[1])
                        end = int(sep[2])
                        # Se coloca la hebra en el formato adecuado
                        if sep[3] == "+":
                            sense = +1
                        else:
                            sense = -1

                    # Se encuentra la correspondencia entre la isla del output y del .gbk
                    if gi == record.id:
                        # Se crea una lista con las coordenadas del CDS
                        spos = sf.ExactPosition(start - 1)
                        epos = sf.ExactPosition(end)
                        l1 = sf.FeatureLocation(spos, epos, strand=sense)

                        # Se lee cada entrada del genbank
                        for rf in record.features:

                            # Asegurarse que el CDS coincida con el del output
                            if l1 == rf.location:
                                # Modificar cada diccionario para agregar inference, product y note
                                pre_inference_list = sep[5]
                                inference_list = list(pre_inference_list.split(','))
                                rf.qualifiers['inference'] = inference_list
                                rf.qualifiers['product'] = sep[6]

                                # Probar si existe el qualifier /note,
                                try:
                                    # Modificar cada diccionario para agregar solamente el /note
                                    list_note = rf.qualifiers['note']
                                    list_note.append("identified in CARD database")
                                    rf.qualifiers['note'] = list_note
                                    c_cds = c_cds + 1
                                except KeyError:
                                    # Si no existe /note, se crea
                                    rf.qualifiers['note'] = "identified in CARD database"
                                    c_cds = c_cds + 1
                            else:
                                continue

            # Escribir la información en el archivo gbk que agrupa toda la información
            with open(out_prokka_gbk, "a") as gbk_output:
                SeqIO.write(record, gbk_output, "genbank")

# Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' % (in_file_short, in_file_short), 'a') as txt_output:
    txt_output.write("   " + str(c_cds) + " CDS annotated with CARD database.\n")

print(str(c_cds) + " CDS annotated with CARD database were included in the data.")

# 2.3.5 Se incorpora en el qualifier /note los CDS anotados por PROKKA para
#       la quinta base de datos que es BacMet2 y aunque cambie la anotación
#       más adelante, quede registro que fue también anotado con BacMet2.

# Contador de CDS incluidos a los datos
c_cds = 0

# En este caso, de la cuarta modificación usaremos como base los archivos CARD
# y le sumaremos la anotación de los resultados BacMet2.
for site in non_redundant_sites:
    # Archivo del cual se obtendrá la información usada como base para la anotación
    in_prokka_gbk = "./%s/results_prokka/%s_GIs/replace_3_from_VFDB_to_CARD.gb" % (in_file_short, site)
    # Archivo con las anotaciones tabuladas que fueron realizadas con BacMet2
    in_data = "./%s/outputs/output_BacMet2.tsv" % (in_file_short)

    # Crear el archivo nuevo, sobreescribiendo si ya existía
    record = ""
    out_prokka_gbk = "./%s/results_prokka/%s_GIs/replace_4_from_CARD_to_BacMet2.gb" % (in_file_short, site)
    with open(out_prokka_gbk, "w") as gbk_output:
        SeqIO.write(record, gbk_output, "genbank")

    # Analizar los archivos de cada carpeta
    with open(in_prokka_gbk, "r") as annot_file:
        # Se le pide a SeqIO que manipule la info del archivo como un genbank
        for record in SeqIO.parse(annot_file, "genbank"):
            # El objeto record es la secuencia y las anotaciones que ya tiene
            # La sequencia está en record.seq y las anotaciones en record.features

            # Abrir la tabla de las anotaciones
            with open(in_data, "r") as data_genes:
                for line in data_genes:
                    sep = line.strip().split("\t")

                    if sep[0] == "tDNA":
                        # Se salta esta línea
                        continue

                    else:
                        # Se procesa la información de esta línea
                        gi = sep[0]
                        start = int(sep[1])
                        end = int(sep[2])
                        # Se coloca la hebra en el formato adecuado
                        if sep[3] == "+":
                            sense = +1
                        else:
                            sense = -1

                    # Se encuentra la correspondencia entre la isla del output y del .gbk
                    if gi == record.id:
                        # Se crea una lista con las coordenadas del CDS
                        spos = sf.ExactPosition(start - 1)
                        epos = sf.ExactPosition(end)
                        l1 = sf.FeatureLocation(spos, epos, strand=sense)

                        # Se lee cada entrada del genbank
                        for rf in record.features:

                            # Asegurarse que el CDS coincida con el del output
                            if l1 == rf.location:
                                # Modificar cada diccionario para agregar inference, product y note
                                pre_inference_list = sep[5]
                                inference_list = list(pre_inference_list.split(','))
                                rf.qualifiers['inference'] = inference_list
                                rf.qualifiers['product'] = sep[6]

                                # Probar si existe el qualifier /note,
                                try:
                                    # Modificar cada diccionario para agregar solamente el /note
                                    list_note = rf.qualifiers['note']
                                    list_note.append("identified in BacMet2 predicted database")
                                    rf.qualifiers['note'] = list_note
                                    c_cds = c_cds + 1
                                except KeyError:
                                    # Si no existe /note, se crea
                                    rf.qualifiers['note'] = "identified in BacMet2 predicted database"
                                    c_cds = c_cds + 1
                            else:
                                continue

            # Escribir la información en el archivo gbk que agrupa toda la información
            with open(out_prokka_gbk, "a") as gbk_output:
                SeqIO.write(record, gbk_output, "genbank")

# Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' % (in_file_short, in_file_short), 'a') as txt_output:
    txt_output.write("   " + str(c_cds) + " CDS annotated with BacMet2 predicted database.\n\n")

print(str(c_cds) + " CDS annotated with BacMet2 predicted database were included in the data.")
print(' ')

# Fin de la segunda parte

actual_time = datetime.datetime.now().strftime('%H:%M:%S')
print('[' + actual_time + '] ======== Paso 2: Anotación general con diferentes bases de datos. FINALIZADO =======')
print(' ')

##############################################################################
########### 3. EJECUTAR PhageBoost (DEBE ESTAR ACTIVO EL AMBIENTE) ###########
##############################################################################

# Inicio de la tercera parte

print(' ')
actual_time = datetime.datetime.now().strftime('%H:%M:%S')
print('[' + actual_time + '] =============== Paso 3: Identificación de probables fagos. INICIANDO ===============')
print(' ')

# 3.1 EJECUCIÓN: PhageBoost toma el multifasta y lo procesa

# No olvidar que se debe ejecutar dentro del ambiente wapi
cmd_phageboost = "PhageBoost -f %s -o ./%s/results_phageboost" % (in_file, in_file_short)
ex_phageboost = subprocess.run(cmd_phageboost, shell=True)

# 3.2 PROCESAMIENTO E INCORPORACIÓN: agregar todos los CDS reconocidos como
#     parte de un fago, según Phageboost, en un archivo genbank para que
#     recopile esto junto a las anotaciones de PROKKA

# Contador de CDS incluidos a los datos
c_cds = 0
c_regions = 0

# Se agregará en cada qualifier /note si el CDS es parte de una región
# compatible con fago. Se analiza el último archivo .gb de las anotaciones de
# PROKKA para agregarle la información correspondiente del archivo output .gff3
for site in non_redundant_sites:

    # Se indica el archivo original desde donde se sumará la información de esta sección
    # ¡¡ATENTO SI CAMBIO EL NOMBRE DEL ÚLTIMO ARCHIVO DE PROKKA!!
    in_phageboost_gbk = "./%s/results_prokka/%s_GIs/replace_4_from_CARD_to_BacMet2.gb" % (in_file_short, site)
    # Archivo de los salida con la predicción de los fagos
    in_data = "./%s/results_phageboost/genecalls_%s.gff3" % (in_file_short, in_file_short)
    in_data_regions = "./%s/results_phageboost/phages_%s.gff" % (in_file_short, in_file_short)

    # Crear el archivo nuevo, sobreescribiendo si ya existía
    record = ""
    out_phageboost_gbk = "./%s/results_prokka/%s_GIs/post_phageboost.gb" % (in_file_short, site)
    with open(out_phageboost_gbk, "w") as gbk_output:
        SeqIO.write(record, gbk_output, "genbank")

    # Analizar los archivos de cada carpeta
    with open(in_phageboost_gbk, "r") as annot_file:
        # Se le pide a SeqIO que manipule la info del archivo como un genbank
        for record in SeqIO.parse(annot_file, "genbank"):
            # El objeto record es la secuencia y las anotaciones que ya tiene
            # La sequencia está en record.seq y las anotaciones en record.features

            # 3.2.1 Abrir la tabla con los resultados de los CDS
            with open(in_data, "r") as data_genes:
                for line in data_genes:
                    sep = line.strip().split("\t")

                    if sep[0] == "contig":
                        # Se salta esta línea
                        continue

                    else:
                        # Se procesa la información de esta línea
                        gi = sep[1]
                        start = int(sep[3])
                        end = int(sep[4])
                        # Se coloca la hebra en el formato adecuado
                        if sep[5] == "1":
                            sense = +1
                        else:
                            sense = -1

                    # Recoger solo los genes de las islas coincidentes y que hayan
                    # sido identificados como fagos.
                    if gi == record.id and "phage" in sep[12]:
                        # Se crea una lista con las coordenadas del CDS
                        spos = sf.ExactPosition(start - 1)
                        epos = sf.ExactPosition(end)
                        l1 = sf.FeatureLocation(spos, epos, strand=sense)

                        # Se lee cada entrada del genbank
                        for rf in record.features:

                            # Asegurarse que el CDS coincida con el del .gff
                            if l1 == rf.location:

                                # Probar si existe el qualifier /note
                                try:
                                    # Modificar cada diccionario para agregar solamente el /note
                                    list_note = rf.qualifiers['note']
                                    list_note.append("identified by PhageBoost as a gene within " + sep[12])
                                    rf.qualifiers['note'] = list_note
                                    c_cds = c_cds + 1
                                except KeyError:
                                    # Si no existe /note, se crea
                                    rf.qualifiers['note'] = "identified by PhageBoost as a gene within " + sep[12]
                                    c_cds = c_cds + 1
                            else:
                                continue
                    else:
                        continue

            # 3.2.2 Abrir la tabla con los resultados de las regiones
            with open(in_data_regions, "r") as data_genes:
                for line in data_genes:
                    sep = line.strip().split("\t")

                    if sep[0] == "##gff-version":
                        # Se salta esta línea
                        continue
                    elif sep[0] == "#seqid":
                        # Se salta esta línea
                        continue

                    else:
                        gi = sep[0]

                        # Solo si el record tiene una región dentro del .gff
                        if gi == record.id:
                            c_regions = c_regions + 1
                            # Se procesa la información de esta línea
                            start = int(sep[3])
                            end = int(sep[4])
                            sense = +1
                            # Se crea una lista con las coordenadas de la región
                            l1 = sf.FeatureLocation(start - 1, end, strand=sense)
                            # Separar la última columna
                            last_column = sep[8]
                            last_column_sp = last_column.split(";")
                            for item in last_column_sp:
                                att = item.split("=")

                                # Eliminar el salto de línea
                                if att[0] == "n_genes":
                                    n_genes = att[1]
                                elif att[0] == "phage_id":
                                    name_phage = att[1]
                                else:
                                    continue
                            f1 = sf.SeqFeature(l1, strand=sense, type="misc_feature", qualifiers={"label": name_phage,
                                                                                                  "note": "This region, predicted by PhageBoost as " + name_phage + ", has " + n_genes + " genes"})
                            record.features.append(f1)

                        else:
                            continue

            # 3.2.3 Escribir la información en el archivo gbk que agrupa toda la información
            with open(out_phageboost_gbk, "a") as gbk_output:
                SeqIO.write(record, gbk_output, "genbank")

# Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' % (in_file_short, in_file_short), 'a') as txt_output:
    txt_output.write("3. PROPHAGES: Phageboost results were " + str(c_cds) + " annotated CDS in " + str(
        c_regions) + " regions. Every region has an independent and named feature.\n   All the information is in the corresponding /note qualifier.\n\n")

print(' ')
print(str(c_cds) + " CDS annotated and " + str(
    c_regions) + " regions, predicted by PhageBoost, were included in the data.")
print(' ')

# Fin de la tercera parte

actual_time = datetime.datetime.now().strftime('%H:%M:%S')
print('[' + actual_time + '] =============== Paso 3: Identificación de probables fagos. FINALIZADO ==============')
print(' ')

##############################################################################
########### 4. EJECUTAR antiSMASH (DEBE ESTAR ACTIVO EL AMBIENTE) ############
##############################################################################

# Inicio de la cuarta parte

print(' ')
actual_time = datetime.datetime.now().strftime('%H:%M:%S')
print('[' + actual_time + '] =========== Paso 4: Identificación de clusters biosintéticos. INICIANDO  ===========')
print(' ')

# 4.1 EJECUCIÓN: antiSMASH  toma el multifasta y lo procesa

# No olvidar que se debe ejecutar dentro del ambiente wapi
cmd_antismash = "antismash %s --output-dir ./%s/results_antismash --cb-knownclusters --genefinding-tool prodigal" % (
in_file, in_file_short)
ex_antismash = subprocess.run(cmd_antismash, shell=True)

# 4.2 PROCESAMIENTO E INCORPORACIÓN: agregar todos los CDS reconocidos como
#     parte de un cluster biosintético, según antiSMASH, en un archivo genbank
#     para que recopile esto junto a las anotaciones de PROKKA y Phageboost. A
#     diferencia de los casos anteriores aquí se debe trabajar con dos genbanks.


# Contador de CDS incluidos a los datos
c_cds = 0

# Primero hay que abrir todos los archivos .gbk, pero como tienen diversos nombres,
# se abrirán, excepto el archivo original de las islas y se agregarán a una lista.
gbkfiles = []
for file in glob.glob("./%s/results_antismash/*.gbk" % in_file_short):
    if in_file_short + ".gbk" in file:
        continue
    else:
        gbkfiles.append(file)

# Primero se abren los .gb de cada una de las islas, luego se revisa cada .gbk
# de salida de antiSMASH para incorporar las regiones correspondientes.
for site in non_redundant_sites:
    # Se indica el archivo original desde donde se sumará la información de esta sección
    in_antismash_gbk = "./%s/results_prokka/%s_GIs/post_phageboost.gb" % (in_file_short, site)

    # Crear el archivo nuevo, sobreescribiendo si ya existía
    record = ""
    out_antismash_gbk = "./%s/results_prokka/%s_GIs/post_antismash.gb" % (in_file_short, site)
    with open(out_antismash_gbk, "w") as gbk_output:
        SeqIO.write(record, gbk_output, "genbank")

    # Analizar los archivos de cada carpeta
    with open(in_antismash_gbk, "r") as annot_file:
        # Se le pide a SeqIO que manipule la info del archivo como un genbank
        for record in SeqIO.parse(annot_file, "genbank"):
            # El objeto record es la secuencia y las anotaciones que ya tiene
            # La sequencia está en record.seq y las anotaciones en record.features

            # Ahora vamos viendo de a uno en uno los outputs de las regiones
            # identificadas por antiSMASH.
            for file in gbkfiles:
                in_data = file

                # Ahora abrimos cada archivo/región de antiSMASH
                with open(in_data, "r") as data_genes:

                    # Por si hubiera más de un record dentro del archivo/región antiSMASH
                    for record_aS in SeqIO.parse(data_genes, "genbank"):
                        # Hacer match del mismo gen en cada archivo
                        if record_aS.description == record.id:
                            for rf_aS in record_aS.features:

                                # Se toma solo la información de los features llamados "proto_core"
                                if rf_aS.type == "proto_core":
                                    # Se toman las coordenadas del archivo procesado por antiSMASH y, como son relativas
                                    # a la región seleccionada, hay que volver a colocar los valores originales
                                    spos = sf.ExactPosition(rf_aS.location.start + int(
                                        record_aS.annotations['structured_comment']['antiSMASH-Data']['Orig. start']))
                                    epos = sf.ExactPosition(rf_aS.location.end + int(
                                        record_aS.annotations['structured_comment']['antiSMASH-Data']['Orig. start']))
                                    sense = rf_aS.location.strand
                                    # Se crea una lista con las coordenadas del CDS
                                    l1 = sf.FeatureLocation(spos, epos, strand=sense)

                                    # Se lee cada entrada del genbank
                                    for rf in record.features:

                                        # Asegurase que coincida con el feature del archivo post_phageboost
                                        if l1 == rf.location:

                                            # Probar si existe el qualifier /note
                                            try:
                                                # Modificar cada diccionario para agregar solamente el /note
                                                list_note = rf.qualifiers['note']
                                                list_note.append("identified by antiSMASH as a " + str(
                                                    rf_aS.qualifiers["product"][0]) + " gene")
                                                rf.qualifiers['note'] = list_note
                                                c_cds = c_cds + 1
                                            except KeyError:
                                                # Si no existe /note, se crea
                                                rf.qualifiers['note'] = "identified by antiSMASH as a " + str(
                                                    rf_aS.qualifiers["product"][0]) + " gene"
                                                c_cds = c_cds + 1
                                        else:
                                            continue

            # Escribir la información en el archivo gbk que agrupa toda la información
            with open(out_antismash_gbk, "a") as gbk_output:
                SeqIO.write(record, gbk_output, "genbank")

# Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' % (in_file_short, in_file_short), 'a') as txt_output:
    txt_output.write("4. BIOSYNTHETIC CLUSTERS: antiSMASH results were " + str(
        c_cds) + " annotated CDS.\n   All the information is in the corresponding /note qualifier.\n\n")

print(str(c_cds) + " CDS annotated with prediction by antiSMASH were included in the data.")
print(' ')

# Fin de la cuarta parte

actual_time = datetime.datetime.now().strftime('%H:%M:%S')
print('[' + actual_time + '] =========== Paso 4: Identificación de clusters biosintéticos. FINALIZADO ===========')
print(' ')

##############################################################################
###### 5. EJECUTAR eggNOG-mapper v2.1.6 (DEBE ESTAR ACTIVO EL AMBIENTE) ######
##############################################################################

# Inicio de la quinta parte

print(' ')
actual_time = datetime.datetime.now().strftime('%H:%M:%S')
print('[' + actual_time + '] ============ Paso 5: Anotación funcional de las secuencias. INICIANDO  =============')
print(' ')

# 5.0 Crear carpeta que albergue todos los resultados de eggNOG

path_eggnog = "results_eggnog"
path = os.path.join(path_file, path_eggnog)
try:
    os.mkdir(path)
except OSError:
    print("Directory /%s not created, because already exists." % path)
else:
    print("Successfully created the directory /%s." % path)
print(' ')

# 5.1 Crear .gbk con todos los CDS calificados como "hypothetical protein" por PROKKA


# Crear el archivo nuevo, sobreescribiendo si ya existía
record = ""
out_pre_eggnog_gbk = "./%s/results_eggnog/hypoprots_%s.gb" % (in_file_short, in_file_short)
with open(out_pre_eggnog_gbk, "w") as gbk_output:
    SeqIO.write(record, gbk_output, "genbank")

# Recorrer todos los archivos post_antiSMASH y quedarse solo con los CDS que
# no hayan sido anotados con proteínas conocidas, es decir, quedarse con los
# digan "hypothetical protein" en el qualifier /product.
for site in non_redundant_sites:
    # Abrir los archivos en cada carpeta
    in_pre_eggnog_gbk = "./%s/results_prokka/%s_GIs/post_antismash.gb" % (in_file_short, site)

    # Analizar los archivos de cada carpeta
    with open(in_pre_eggnog_gbk, "r") as annot_file:
        # Se le pide a SeqIO que manipule la info del archivo como un genbank
        for record in SeqIO.parse(annot_file, "genbank"):
            # El objeto record es la secuencia y las anotaciones que ya tiene
            # La sequencia está en record.seq y las anotaciones en record.features

            # Crear un record con la misma información para ir agregando lo que interesa
            record_pre_eggnog = SeqRecord(record.seq)
            record_pre_eggnog.id = record.id
            record_pre_eggnog.description = record.description
            record_pre_eggnog.name = record.name
            record_pre_eggnog.annotations = record.annotations
            record_pre_eggnog.features = []

            # Se lee cada entrada del genbank
            for rf in record.features:
                if rf.type == "source":
                    # Agregarlo a record_pre_eggnog
                    record_pre_eggnog.features.append(rf)
                elif rf.type == "CDS":
                    if rf.qualifiers['product'][0] == "hypothetical protein":
                        # Agregar a un record nuevo que será escrito en out_pre_eggnog_gbk
                        record_pre_eggnog.features.append(rf)
                    else:
                        continue
                else:
                    continue

            # Escribir la información del record_pre_eggnog en el archivo post_eggnog
            with open(out_pre_eggnog_gbk, "a") as gbk_output:
                SeqIO.write(record_pre_eggnog, gbk_output, "genbank")

# 5.2 Transformar el .gbk en .fasta de aminoácidos para que lo pueda procesar eggNOG

# Crear el archivo nuevo, sobreescribiendo si ya existía
record = ""
out_pre_eggnog_fasta = "./%s/results_eggnog/hypoprots_%s.faa" % (in_file_short, in_file_short)
with open(out_pre_eggnog_fasta, "w") as output_pre_eggnog:
    SeqIO.write(record, output_pre_eggnog, "fasta")

# Abrir el archivo con la información en .gbk
output_handle = open(out_pre_eggnog_fasta, "w")
c_eggnog = 0
with open(out_pre_eggnog_gbk, "r") as input_pre_eggnog:
    for record in SeqIO.parse(input_pre_eggnog, "genbank"):
        for seq_feature in record.features:

            # Abrir y escribir el archivo fasta de aminoácidos
            try:
                output_handle.write(">%s;%s;[%s:%s](%s)\n%s\n" % (
                seq_feature.qualifiers['locus_tag'][0], record.name, int(seq_feature.location.start) + 1,
                seq_feature.location.end, seq_feature.location.strand, seq_feature.qualifiers['translation'][0]))
                c_eggnog = c_eggnog + 1
            except KeyError:
                continue
output_handle.close()

# Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' % (in_file_short, in_file_short), 'a') as txt_output:
    txt_output.write(
        "5. FUNCTIONAL ANNOTATION: eggNOG used a reduced file containing only 'hypothetical proteins'.\n   The file had " + str(
            c_eggnog) + " features.\n")

print(str(c_eggnog) + " hypothetical proteins were included in the .faa file.")
print(' ')

# 5.3 EJECUCIÓN: Tomar el .fasta de aminoácidos con la colección reducida

cmd_eggnog = "python ~/anaconda3/envs/wapi/bin/emapper.py -i ./%s/results_eggnog/hypoprots_%s.faa --output_dir ./%s/results_eggnog -o hypoprots_%s --cpu 0" % (
in_file_short, in_file_short, in_file_short, in_file_short)
ex_eggnog = subprocess.run(cmd_eggnog, shell=True)

# 5.4 PROCESAMIENTO: Crear el output con los datos relevantes de eggNOG

# Crear el archivo nuevo, sobreescribiendo si ya existía
header_eggnog_tsv = ["Island", "Locus Tag", "Start", "End", "Strand", "Code|COG Category", "Description",
                     "Preferred Name"]
with open("./%s/outputs/output_eggnog.tsv" % in_file_short, 'w', newline='') as f_output:
    tsv_output = csv.writer(f_output, delimiter='\t')
    tsv_output.writerow(header_eggnog_tsv)

    # 5.4.1 Importar el diccionario de categorías.
with open("dict_eggnog.txt") as f_dict:
    data_dict_COG = f_dict.read()
dict_COG_cat = json.loads(data_dict_COG)

# Tomar el archivo de las anotaciones hechas por eggNOG
eggnog_in_file = "./%s/results_eggnog/hypoprots_%s.emapper.annotations" % (in_file_short, in_file_short)
COG_categories = []
list_of_genes = []

# Analizar el archivo de output de los resultados de eggNOG
with open(eggnog_in_file, 'r') as data_eggnog:
    for line in data_eggnog:

        if '#' in line:
            # Se salta esta línea
            continue

        else:
            # 5.4.2 Se procesa la información de esta línea
            gene_eggnog = line.split("\t")
            query_column = gene_eggnog[0]
            query_column_sep = query_column.split(";")
            # Tomar el nombre de la isla y el locus tag
            island_column = query_column_sep[1]
            locus_tag_column = query_column_sep[0]

            # 5.4.3 Tomar las coordenadas y crear variables para poder separar la información
            coord_column = query_column_sep[2]
            pos_sep_numbers = coord_column.rfind(":")
            pos_sep_strand = coord_column.rfind("]")
            # Dejar la info de las coordenadas en 3 columnas diferentes
            start_column = coord_column[1:pos_sep_numbers]
            end_column = coord_column[pos_sep_numbers + 1:pos_sep_strand]
            strand_column = coord_column[pos_sep_strand + 2:-1]

            # 5.4.4 Crear la columna con el nombre de las categorías.
            COG_column = gene_eggnog[6]
            # Lista especial que contendrá solo un caracter por elemento, usada más adelante.
            COG_column_sep = COG_column


            def split(COG_column_sep):
                return [char for char in COG_column_sep]


            COG_categories.extend(COG_column_sep)

            # 5.4.5 Tomar las categorías COG, separar las que tengan más de una categoría (1, 2, 3 o más).
            if len(COG_column) == 1:
                category_name_column = dict_COG_cat[COG_column[0]]
            elif len(COG_column) == 2:
                category_name_column = dict_COG_cat[COG_column[0]] + "; " + dict_COG_cat[COG_column[1]]
            elif len(COG_column) == 3:
                category_name_column = dict_COG_cat[COG_column[0]] + "; " + dict_COG_cat[COG_column[1]] + "; " + \
                                       dict_COG_cat[COG_column[2]]
            elif len(COG_column) > 3:
                category_name_column = dict_COG_cat[COG_column[0]] + "; " + dict_COG_cat[COG_column[1]] + "; " + \
                                       dict_COG_cat[COG_column[2]] + "; and more categories..."
                print('"The query " query_column " has 4 or more categories, will be shown only the three first"')
            else:
                continue

            # 5.4.6 Tomar la columna de la descripción y del nombre preferido
            description_column = gene_eggnog[7]
            preferred_name_column = gene_eggnog[8]

            # 5.4.7 Se junta toda la información útil en una variable
            useful_info = [island_column, locus_tag_column, start_column, end_column, strand_column,
                           COG_column + "|" + category_name_column, description_column, preferred_name_column]
            list_of_genes.append(useful_info)

            # 5.4.8 Escribir la información en el archivo output
with open("./%s/outputs/output_eggnog.tsv" % in_file_short, 'a+', newline='') as f_output:
    tsv_output = csv.writer(f_output, delimiter='\t')
    tsv_output.writerows(list_of_genes)

# 5.5 INCORPORACIÓN: agregar toda la información del output en un archivo .gb
#     para que recopile toda la información obtenida hasta el momento.

# Varios contadores
c_cds = 0
c_hypo_prots = 0
c_rec = 0

# Se leerá cada archivo de salida de antismash para colocar en /note si el CDS
# de "hypothetical protein" tiene una anotación otorgada por eggNOG
for site in non_redundant_sites:
    # Archivo del cual se obtendrá la información usada como base para la anotación
    in_eggnog_gbk = "./%s/results_prokka/%s_GIs/post_antismash.gb" % (in_file_short, site)
    # Archivo con las anotaciones tabuladas que fueron realizadas con eggNOG
    in_data = "./%s/outputs/output_eggnog.tsv" % (in_file_short)

    # Crear el archivo nuevo, sobreescribiendo si ya existía
    record = ""
    out_eggnog_gbk = "./%s/results_prokka/%s_GIs/post_eggnog.gb" % (in_file_short, site)
    with open(out_eggnog_gbk, "w") as gbk_output:
        SeqIO.write(record, gbk_output, "genbank")

    # Analizar los archivos de cada carpeta
    with open(in_eggnog_gbk, "r") as annot_file:
        # Se le pide a SeqIO que manipule la info del archivo como un genbank
        for record in SeqIO.parse(annot_file, "genbank"):
            # Contar el numero de records
            c_rec = c_rec + 1
            # El objeto record es la secuencia y las anotaciones que ya tiene
            # La sequencia está en record.seq y las anotaciones en record.features

            # Abrir la tabla de las anotaciones
            with open(in_data, "r") as data_genes:
                for line in data_genes:
                    sep = line.strip().split("\t")

                    if sep[0] == "Island":
                        # Se salta esta línea
                        continue

                    elif sep[5] == "-|Not recognized by eggNOG-mapper":
                        # Como no tiene información relevante, se salta esta línea
                        c_hypo_prots = c_hypo_prots + 1
                        continue

                    else:
                        # Se procesa la información de esta línea
                        gi = sep[0]
                        start = int(sep[2])
                        end = int(sep[3])
                        # Se coloca la hebra en el formato adecuado
                        if sep[4] == "1":
                            sense = +1
                        else:
                            sense = -1
                        # Se modifica el nombre de la categoría para un mejor orden
                        sep[5] = sep[5].replace("|", " (")

                    # Se encuentra la correspondencia entre la isla del output y del .tsv
                    if gi == record.id:
                        # Se crea una lista con las coordenadas del CDS
                        spos = sf.ExactPosition(start - 1)
                        epos = sf.ExactPosition(end)
                        l1 = sf.FeatureLocation(spos, epos, strand=sense)

                        # Se lee cada entrada del genbank
                        for rf in record.features:

                            # Asegurarse que el CDS coincida con el del output
                            if l1 == rf.location:

                                # Dependiendo de si tiene nombre preferido, se escriben dos cosas diferentes en /note
                                if sep[7] == "-":
                                    # Probar si existe el qualifier /note
                                    try:
                                        # Modificar cada diccionario para agregar solamente el /note
                                        list_note = rf.qualifiers['note']
                                        list_note.append("identified by eggNOG as a gene in the category " + sep[
                                            5] + "). Its description is '" + sep[6] + "'")
                                        rf.qualifiers['note'] = list_note
                                        c_cds = c_cds + 1
                                    except KeyError:
                                        # Si no existe /note, se crea
                                        rf.qualifiers['note'] = "identified by eggNOG as a gene in the category " + sep[
                                            5] + "). Its description is '" + sep[6] + "'"
                                        c_cds = c_cds + 1

                                else:
                                    # Probar si existe el qualifier /note
                                    try:
                                        # Modificar cada diccionario para agregar solamente el /note
                                        list_note = rf.qualifiers['note']
                                        list_note.append(
                                            "identified by eggNOG as the '" + sep[7] + "' gene in the category " + sep[
                                                5] + "). Its description is '" + sep[6] + "'")
                                        rf.qualifiers['note'] = list_note
                                        c_cds = c_cds + 1
                                    except KeyError:
                                        # Si no existe /note, se crea
                                        rf.qualifiers['note'] = "identified by eggNOG as the '" + sep[
                                            7] + "' gene in the category " + sep[5] + "). Its description is '" + sep[
                                                                    6] + "'"
                                        c_cds = c_cds + 1
                            else:
                                continue

            # Escribir la información en el archivo gbk que agrupa toda la información
            with open(out_eggnog_gbk, "a") as gbk_output:
                SeqIO.write(record, gbk_output, "genbank")

# Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' % (in_file_short, in_file_short), 'a') as txt_output:
    txt_output.write("   Of these, " + str(c_cds) + " 'hypothetical proteins' were annotated and " + str(
        int(c_hypo_prots / c_rec)) + " CDS could not be annotated.\n")

print(str(c_cds) + " 'hypothetical proteins' annotated by eggNOG were included in the data.")
print("Unfortunately, " + str(
    int(c_hypo_prots / c_rec)) + " CDS could not be annotated by a database used in PROKKA or another program used.")
print(' ')

# 5.6 Como es interesante analizar las integrasas encontradas por eggNOG, se
# creará un archivo gbk especial que contendrá solo estos CDS.

# Contador de islas genómicas incluidas en el archivo final (c de count)
c_integrases = 0

# Crear el archivo nuevo, sobreescribiendo si ya existía
record = ""
out_integrase_gbk = "./%s/%s_integrases.gbk" % (in_file_short, in_file_short)
with open(out_integrase_gbk, "w") as gbk_output:
    SeqIO.write(record, gbk_output, "genbank")

# Se recorren todas las carpetas
for site in non_redundant_sites:
    in_integrase_gbk = "./%s/results_prokka/%s_GIs/post_eggnog.gb" % (in_file_short, site)

    # Analizar los archivos de cada carpeta
    with open(in_integrase_gbk, "r") as annot_file:
        # Se le pide a SeqIO que manipule la info del archivo como un genbank
        for record in SeqIO.parse(annot_file, "genbank"):
            # El objeto record es la secuencia y las anotaciones que ya tiene
            # La sequencia está en record.seq y las anotaciones en record.features

            # Aquí ingresarán solo los features que tengan un CDS de integrasa
            # que haya sido reconocido con eggNOG
            record_integrase = []

            # Crear un record temporal con la misma información del archivo post_eggnog
            record_putative = SeqRecord(record.seq)
            record_putative.id = record.id
            record_putative.description = record.description
            record_putative.name = record.name
            record_putative.annotations = record.annotations
            record_putative.features = []

            # Se lee cada entrada del genbank
            for rf in record.features:
                if rf.type == "source":
                    # Agregarlo a record_putative mientras se espera si reconoce algún CDS
                    record_putative.features.append(rf)

                # Hay muchos features que no tienen el qualifier /note, primero,
                # hay que probar. Luego si "integrase" y "eggNOG" están en el mismo /note
                try:
                    # Buscar las dos palabras en cada elemento de la lista ""note
                    if "integrase" in rf.qualifiers["note"][0]:
                        if "eggNOG" in rf.qualifiers["note"][0]:
                            # Agregar el feature al record temporal
                            record_putative.features.append(rf)
                            c_integrases = c_integrases + 1
                            # Se envía al record que se escribirá en el archivo
                            record_integrase = record_putative
                    # Aunque debiese ser el primer elemento en /note,
                    # igual buscará en otros dos, por si ya tuviera de antes otras anotaciones
                    elif "integrase" in rf.qualifiers["note"][1]:
                        if "eggNOG" in rf.qualifiers["note"][1]:
                            # Agregar el feature al record temporal
                            record_putative.features.append(rf)
                            c_integrases = c_integrases + 1
                            # Se envía al record que se escribirá en el archivo
                            record_integrase = record_putative
                    elif "integrase" in rf.qualifiers["note"][2]:
                        if "eggNOG" in rf.qualifiers["note"][2]:
                            # Agregar el feature al record temporal
                            record_putative.features.append(rf)
                            c_integrases = c_integrases + 1
                            # Se envía al record que se escribirá en el archivo
                            record_integrase = record_putative
                    elif "integrase" in rf.qualifiers["note"][3]:
                        if "eggNOG" in rf.qualifiers["note"][3]:
                            # Agregar el feature al record temporal
                            record_putative.features.append(rf)
                            c_integrases = c_integrases + 1
                            # Se envía al record que se escribirá en el archivo
                            record_integrase = record_putative
                    elif "integrase" in rf.qualifiers["note"][4]:
                        if "eggNOG" in rf.qualifiers["note"][4]:
                            # Agregar el feature al record temporal
                            record_putative.features.append(rf)
                            c_integrases = c_integrases + 1
                            # Se envía al record que se escribirá en el archivo
                            record_integrase = record_putative
                    elif "integrase" in rf.qualifiers["note"][5]:
                        if "eggNOG" in rf.qualifiers["note"][5]:
                            # Agregar el feature al record temporal
                            record_putative.features.append(rf)
                            c_integrases = c_integrases + 1
                            # Se envía al record que se escribirá en el archivo
                            record_integrase = record_putative
                except IndexError:
                    # Si no existe ese /note , se ingnora
                    continue
                except KeyError:
                    # Si no existe /note, se ingnora
                    continue

            # Escribir la información del record_integrase en el archivo de las integrasas
            with open(out_integrase_gbk, "a") as gbk_output:
                SeqIO.write(record_integrase, gbk_output, "genbank")

# Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' % (in_file_short, in_file_short), 'a') as txt_output:
    txt_output.write(
        "   Furthermore, the information of annotation made it possible to create a file containing " + str(
            c_integrases) + " integrases predicted by eggNOG.\n\n")

print(str(c_integrases) + " integrases predicted in eggNOG were included in %s_integrases.gbk." % (in_file_short))
print(' ')

# Fin de la quinta parte

actual_time = datetime.datetime.now().strftime('%H:%M:%S')
print('[' + actual_time + '] ============ Paso 5: Anotación funcional de las secuencias. FINALIZADO =============')
print(' ')

##############################################################################
############## 6. BUSCAR SECUENCIAS CORTAS DE DNA EN TADB Y oriT #############
##############################################################################

# Inicio de la sexta parte

print(' ')
actual_time = datetime.datetime.now().strftime('%H:%M:%S')
print('[' + actual_time + '] ========= Paso 6: Búsqueda de secuencias de DNA en TADB y oriT. INICIANDO  =========')
print(' ')

# Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' % (in_file_short, in_file_short), 'a') as txt_output:
    txt_output.write("6. NUCLEOTIDE SEQUENCES: nhmmer used several databases to identify relevant sequences:\n")

# 6.1 Crear la carpeta donde estarán los resultados

path_nhmmer = "results_nhmmer"
path = os.path.join(path_file, path_nhmmer)
try:
    os.mkdir(path)
except OSError:
    print("Directory /%s not created, because already exists." % path)
else:
    print("Successfully created the directory /%s." % path)
print(' ')

# 6.2 EJECUCIÓN: se usan las dos bases de datos por separado en nhmmer

# Análisis con la base de datos de oriT
print("Processing the file with oriTDB database...")
cmd_nhmmer_1 = "nhmmer -o ./%s/results_nhmmer/%s_oriT_output.txt --tblout ./%s/results_nhmmer/%s_oriT_tbl.txt --dfamtblout ./%s/results_nhmmer/%s_oriT_dfam.txt --cpu 100 --dna ./%s ./oriT_all.fas" % (
in_file_short, in_file_short, in_file_short, in_file_short, in_file_short, in_file_short, in_file)
ex_nhmmer_1 = subprocess.run(cmd_nhmmer_1, shell=True)
print("Finished\n")

# Análisis con la base de datos de RNA de TADB
print("Processing the file with TADB's RNAs database...")
cmd_nhmmer_2 = "nhmmer -o ./%s/results_nhmmer/%s_RNA_output.txt --tblout ./%s/results_nhmmer/%s_RNA_tbl.txt --dfamtblout ./%s/results_nhmmer/%s_RNA_dfam.txt --cpu 100 --dna ./%s ./db_TADB_RNA.fasta" % (
in_file_short, in_file_short, in_file_short, in_file_short, in_file_short, in_file_short, in_file)
ex_nhmmer_2 = subprocess.run(cmd_nhmmer_2, shell=True)
print("Finished\n")

# 6.3 PROCESAMIENTO: se toma la información del archivo "tabulado" y se crea
#     el archivo de output correspondiente.

# 6.3.1 Análisis de oriT
# Nombre del archivo de salida de nhmmer
in_data = "./%s/results_nhmmer/%s_oriT_tbl.txt" % (in_file_short, in_file_short)

# Crear el archivo nuevo, con solo el encabezado
header_nhmmer_tsv = ["target name", "accession", "query name", "accession", "hmmfrom", "hmm to", "alifrom", "ali to",
                     "envfrom", "env to", "sq len", "strand", "E-value", "score", "bias", "description of target"]
with open("./%s/outputs/output_oriT.tsv" % (in_file_short), 'w', newline='') as f_output:
    tsv_output = csv.writer(f_output, delimiter='\t')
    tsv_output.writerow(header_nhmmer_tsv)

# Lista que contendrá toda la información que será escrita
all_sep = []

# Abrir el archivo de salida de nhmmer
with open(in_data, "r") as data_genes:
    for line in data_genes:

        if "#" in line:
            continue

        else:
            # Se procesará la línea
            sep = []
            # Se separan las líneas donde haya un doble espacio
            sep = line.strip().split("  ")
            # Se juntarán las líneas con un espacio entremedio,
            # esto eliminará los elementos que son solo un espacio
            sep = " ".join(sep).split()
            last = int(len(sep))
            # Loop que juntará todas las palabras en el elemento de index = 15
            c_word = 15
            list_words_description = []

            while (c_word < int(last)):
                list_words_description.append(sep[c_word])
                new_description = " ".join(list_words_description)
                c_word = c_word + 1

            sep[15] = new_description
            # Se borrarán los elementos sobrantes
            del (sep[16:last + 1])
            all_sep.append(sep)

    # Escribir la información en el archivo nuevo
    with open("./%s/outputs/output_oriT.tsv" % (in_file_short), 'a', newline='') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        tsv_output.writerows(all_sep)

# 6.3.2 Análisis de los RNA de TADB
# Nombre del archivo de salida de nhmmer
in_data = "./%s/results_nhmmer/%s_RNA_tbl.txt" % (in_file_short, in_file_short)

# Crear el archivo nuevo, con solo el encabezado
header_nhmmer_tsv = ["target name", "accession", "query name", "accession", "hmmfrom", "hmm to", "alifrom", "ali to",
                     "envfrom", "env to", "sq len", "strand", "E-value", "score", "bias", "description of target"]
with open("./%s/outputs/output_RNA.tsv" % (in_file_short), 'w', newline='') as f_output:
    tsv_output = csv.writer(f_output, delimiter='\t')
    tsv_output.writerow(header_nhmmer_tsv)

# Lista que contendrá toda la información que será escrita
all_sep = []

# Abrir el archivo de salida de nhmmer
with open(in_data, "r") as data_genes:
    for line in data_genes:

        if "#" in line:
            continue

        else:
            # Se procesará la línea
            sep = []
            # Se separan las líneas donde haya un doble espacio
            sep = line.strip().split("  ")
            # Se juntarán las líneas con un espacio entremedio,
            # esto eliminará los elementos que son solo un espacio
            sep = " ".join(sep).split()
            last = int(len(sep))
            # Loop que juntará todas las palabras en el elemento de index = 15
            c_word = 15
            list_words_description = []

            while (c_word < int(last)):
                list_words_description.append(sep[c_word])
                new_description = " ".join(list_words_description)
                c_word = c_word + 1

            sep[15] = new_description
            # Se borrarán los elementos sobrantes
            del (sep[16:last + 1])
            all_sep.append(sep)

    # Escribir la información en el archivo nuevo
    with open("./%s/outputs/output_RNA.tsv" % (in_file_short), 'a', newline='') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        tsv_output.writerows(all_sep)

# 6.4 INCORPORACIÓN: Tomar la información de los archivos de output de oriT y
#     los RNA de TADB para incorporarlos a un .gbk que acumule toda la
#     información hasta este punto.


# 6.4.1 Información de las secuencias de transferencia oriT

# Contador de CDS incluidos a los datos
c_oriT = 0

# Se creará un feature nuevo para cada secuencia identificada, esto se
# escribirá en el .gbk recopilador de la información
for site in non_redundant_sites:

    # Se indica el archivo original desde donde se sumará la información de esta sección
    in_oriT_gbk = "./%s/results_prokka/%s_GIs/post_eggnog.gb" % (in_file_short, site)
    # Archivo de los salida con la predicción de los fagos
    in_data = "./%s/outputs/output_oriT.tsv" % (in_file_short)

    # Crear el archivo nuevo, sobreescribiendo si ya existía
    record = ""
    out_oriT_gbk = "./%s/results_prokka/%s_GIs/post_oriT.gb" % (in_file_short, site)
    with open(out_oriT_gbk, "w") as gbk_output:
        SeqIO.write(record, gbk_output, "genbank")

    # Analizar los archivos de cada carpeta
    with open(in_oriT_gbk, "r") as annot_file:
        # Se le pide a SeqIO que manipule la info del archivo como un genbank
        for record in SeqIO.parse(annot_file, "genbank"):
            # El objeto record es la secuencia y las anotaciones que ya tiene
            # La sequencia está en record.seq y las anotaciones en record.features

            # Abrir la tabla con los resultados de las secuencias
            with open(in_data, "r") as data_genes:
                for line in data_genes:
                    sep = line.strip().split("\t")

                    if sep[0] == "target name":
                        # Se salta esta línea
                        continue

                    else:
                        gi = sep[2]

                        # Solo si el record tiene una región dentro del .tsv
                        if gi == record.id:
                            c_oriT = c_oriT + 1
                            # Se procesa la información de esta línea
                            start = int(sep[4])
                            end = int(sep[5])
                            if start > end:
                                sense = -1
                            else:
                                sense = +1
                            # Se crea una lista con las coordenadas de la región
                            l1 = sf.FeatureLocation(start - 1, end, strand=sense)

                            # Tomar toda la información del encabezado
                            if "-" == sep[15]:
                                oriT_name = sep[0]
                            else:
                                oriT_name = sep[0] + " " + sep[15]

                            f1 = sf.SeqFeature(l1, strand=sense, type="oriT", qualifiers={"gene": sep[0],
                                                                                          "note": "This oriT sequence was predicted by nhmmer using oriTDB as '" + oriT_name + "', its positions in target sequence are: (" +
                                                                                                  sep[6] + ".." + sep[
                                                                                                      7] + "). For further information, please consult oriTDB webpage"})
                            record.features.append(f1)

                        else:
                            continue

            # Escribir la información en el archivo gbk que agrupa toda la información
            with open(out_oriT_gbk, "a") as gbk_output:
                SeqIO.write(record, gbk_output, "genbank")

# Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' % (in_file_short, in_file_short), 'a') as txt_output:
    txt_output.write("   " + str(c_oriT) + " oriT sequences added with oriTDB database.\n")

print(str(c_oriT) + " oriT sequences with oriTDB database were included in the data.")

# 6.4.2 Información de las antitoxinas RNA de la base de datos TADB

# Contador de CDS incluidos a los datos
c_RNA = 0

# Se creará un feature nuevo para cada secuencia identificada, esto se
# escribirá en el .gbk recopilador de la información
for site in non_redundant_sites:

    # Se indica el archivo original desde donde se sumará la información de esta sección
    in_RNA_gbk = "./%s/results_prokka/%s_GIs/post_oriT.gb" % (in_file_short, site)
    # Archivo de los salida con la predicción de los fagos
    in_data = "./%s/outputs/output_RNA.tsv" % (in_file_short)

    # Crear el archivo nuevo, sobreescribiendo si ya existía
    record = ""
    out_RNA_gbk = "./%s/results_prokka/%s_GIs/post_RNA.gb" % (in_file_short, site)
    with open(out_RNA_gbk, "w") as gbk_output:
        SeqIO.write(record, gbk_output, "genbank")

    # Analizar los archivos de cada carpeta
    with open(in_RNA_gbk, "r") as annot_file:
        # Se le pide a SeqIO que manipule la info del archivo como un genbank
        for record in SeqIO.parse(annot_file, "genbank"):
            # El objeto record es la secuencia y las anotaciones que ya tiene
            # La sequencia está en record.seq y las anotaciones en record.features

            # Abrir la tabla con los resultados de las secuencias
            with open(in_data, "r") as data_genes:
                for line in data_genes:
                    sep = line.strip().split("\t")

                    if sep[0] == "target name":
                        # Se salta esta línea
                        continue

                    else:
                        gi = sep[2]

                        # Solo si el record tiene una región dentro del .tsv
                        if gi == record.id:
                            c_RNA = c_RNA + 1
                            # Se procesa la información de esta línea
                            start = int(sep[4])
                            end = int(sep[5])
                            if start > end:
                                sense = -1
                            else:
                                sense = +1
                            # Se crea una lista con las coordenadas de la región
                            l1 = sf.FeatureLocation(start - 1, end, strand=sense)
                            # Tomar toda la información del encabezado
                            RNA_name = sep[0] + " " + sep[15]
                            f1 = sf.SeqFeature(l1, strand=sense, type="ncRNA",
                                               qualifiers={"gene": "antitoxin RNA", "ncRNA_class": "antisense_RNA",
                                                           "note": "This antisense RNA sequence was predicted by nhmmer using TADB as '" + RNA_name + "', its positions in target sequence are: (" +
                                                                   sep[6] + ".." + sep[7] + ")"})
                            record.features.append(f1)

                        else:
                            continue

            # Escribir la información en el archivo gbk que agrupa toda la información
            with open(out_RNA_gbk, "a") as gbk_output:
                SeqIO.write(record, gbk_output, "genbank")

# Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' % (in_file_short, in_file_short), 'a') as txt_output:
    txt_output.write("   " + str(c_RNA) + " antitoxin RNAs sequences added with TADB database.\n\n")

print(str(c_RNA) + " antitoxin RNAs sequences with TADB database were included in the data.")
print(' ')

# Fin de la sexta parte

actual_time = datetime.datetime.now().strftime('%H:%M:%S')
print('[' + actual_time + '] ========= Paso 6: Búsqueda de secuencias de DNA en TADB y oriT. FINALIZADO  ========')
print(' ')

##############################################################################
################ 7. GENERAR UN ARCHIVO FINAL QUE RESUMA #################
##############################################################################

# Inicio de la séptima parte

print(' ')
actual_time = datetime.datetime.now().strftime('%H:%M:%S')
print('[' + actual_time + '] =========== Paso 7: Generación de un único archivo de salida. INICIANDO  ===========')
print(' ')

# JUNTAR LOS ARCHIVOS .GBK DE CADA ISLA EN UN ÚNICO GENBANK QUE CONTENGA
# TODAS LAS ISLAS NO REDUNDANTES ANALIZADAS. AL MISMO TIEMPO, CREAR ARCHIVOS
# SEPARADOS PARA CADA ISLA CON LA ANOTACIÓN FINAL.

# 7.1 Crear la nueva carpeta de las islas por separado
path_file = in_file_short
path_GIs = "GI_final_annotation_gbk"
path = os.path.join(path_file, path_GIs)
try:
    os.mkdir(path)
except OSError:
    print("Directory /%s not created, because already exists." % path)
else:
    print("Successfully created the directory /%s." % path)
print(' ')

# 7.2 Escribir todos los archivos de salida
# Contador de islas genómicas incluidas en el archivo final (c de count)
c_records = 0

# Crear el archivo nuevo, sobreescribiendo si ya existía
record = ""
out_final_gbk = "./%s/%s_final_annotation.gbk" % (in_file_short, in_file_short)
with open(out_final_gbk, "w") as gbk_output:
    SeqIO.write(record, gbk_output, "genbank")

# Se recorren todas las carpetas
for site in non_redundant_sites:
    in_final_gbk = "./%s/results_prokka/%s_GIs/post_RNA.gb" % (in_file_short, site)

    # Analizar los archivos de cada carpeta
    with open(in_final_gbk, "r") as annot_file:
        # Se le pide a SeqIO que manipule la info del archivo como un genbank
        for record in SeqIO.parse(annot_file, "genbank"):
            # Contar el número de islas incorporadas
            c_records = c_records + 1

            # 7.2.1 Escribir la información de este archivo en el archivo final único
            with open(out_final_gbk, "a") as gbk_output:
                SeqIO.write(record, gbk_output, "genbank")

            # 7.2.2 También escribir la información en archivos individuales
            with open("./%s/GI_final_annotation_gbk/%s_final_annotation.gbk" % (in_file_short, record.id),
                      "w") as gbk_GI_output:
                SeqIO.write(record, gbk_GI_output, "genbank")

# Escribir línea en archivo .txt para el resumen del procesamiento del archivo
with open('./%s/%s_final_annotation.txt' % (in_file_short, in_file_short), 'a') as txt_output:
    txt_output.write("Finally, all the information was added to a unique file: " + str(
        c_records) + " genomic islands were included in %s_final_annotation.gbk. Also, " % (in_file_short) + str(
        c_records) + " files were created for each genomic island in /" + str(path_GIs) + " directory.")

print(str(c_records) + " genomic islands were included in %s_final_annotation.gbk and in single-GI files in /" % (
    in_file_short) + str(path_GIs) + " directory.")
print(' ')

# Fin de la séptima parte

actual_time = datetime.datetime.now().strftime('%H:%M:%S')
print('[' + actual_time + '] =========== Paso 7: Generación de un único archivo de salida. FINALIZADO  ==========')
print(' ')

##############################################################################
############################### FIN DEL SCRIPT ###############################
##############################################################################

# Fin del script

print(' ')
actual_time = datetime.datetime.now().strftime('%H:%M:%S')
print('[' + actual_time + '] ================================= SCRIPT TERMINADO =================================')
print(' ')
print('Gracias por usar este humilde programa.')
print(' ')