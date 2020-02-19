import os
import csv
from itertools import groupby
from biothings.utils.dataload import unlist

def load_data(data_folder):
    gene2disease_file = os.path.join(data_folder, "MGI_DO.rpt")
    gene2phenotype_file = os.path.join(data_folder, "MGI_Geno_DiseaseDO.rpt")
    gene2homolog_file = os.path.join(data_folder, "MGI_Gene_Model_Coord.rpt")
    final_result = {}
    disease_result = load_gene2disese(gene2disease_file)
    phenotype_result = load_gene2phenotype(gene2phenotype_file)
    homolog_result = load_gene2homolog(gene2homolog_file)
    for mgi, info in disease_result.items():
        if mgi not in final_result:
            final_result[mgi] = {"_id": mgi, "mgi": {}}
        final_result[mgi]["mgi"]["associated_with_disease"] = info

    for mgi, phenotype_info in phenotype_result.items():
        if mgi not in final_result:
            final_result[mgi] = {"_id": mgi, "mgi": {}}
        final_result[mgi]["mgi"]["associated_with_phenotype"] = phenotype_info

    for mgi, entrez_info in homolog_result.items():
        if mgi not in final_result:
            final_result[mgi] = {"_id": mgi, "mgi": {}}
        final_result[mgi]["mgi"]["has_homolog"] = entrez_info
    for _doc in final_result.values():
        yield _doc
    

def load_gene2disese(file_path):
    with open(file_path) as f:
        next(f)
        data = list(csv.reader(f, delimiter='\t'))
        result = {}
        data = sorted(data, key=lambda x: x[-1])
        for k, g in groupby(data, lambda x: x[-1]):
            if k:
                res = []
                g = list(g)
                for _item in g:
                    _res = {'doid': _item[0],
                            'name': _item[1],
                            'omim': _item[2].split(':')[-1] if _item[2] else _item[2]}
                    _res = unlist(_res)
                    res.append(_res)
                result[k] = res
    return result


def load_gene2phenotype(file_path):
    with open(file_path) as f:
        next(f)
        phenotype_data = list(csv.reader(f, delimiter='\t'))
        phenotype_data = sorted(phenotype_data, key=lambda x: x[6])
        phenotype_result = {}
        for k, g in groupby(phenotype_data, lambda x: x[6]):
            if k:
                res = []
                g = list(g)
                for _item in g:
                    _res = {'allelic_composition': _item[0],
                            'allele_symbol': _item[1],
                            'allele_id': _item[2],
                            'genetic_background': _item[3],
                            'mp': _item[4],
                            'pubmed': _item[5].split('|')}
                    _res = unlist(_res)
                    res.append(_res)
                phenotype_result[k] = res
    return phenotype_result


def load_gene2homolog(file_path):
    with open(file_path) as f:
        next(f)
        entrez_data = list(csv.reader(f, delimiter='\t'))
        entrez_data = sorted(entrez_data, key=lambda x: x[0])
        entrez_result = {}
        for k, g in groupby(entrez_data, lambda x: x[0]):
            if k:
                entrez_result[k] = [_item[5] for _item in g]
    return entrez_result
