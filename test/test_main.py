import os
import yaml
from cgi_bin import Epitope
from cgi_bin import PHAGE

def load_hla_data(yaml_file):
    lines = []
    with open(yaml_file) as f:
        data = yaml.load(f, Loader=yaml.SafeLoader)
    for hla in data:
        for codon in data[hla]:
            to_add = [hla, codon, data[hla][codon]]
            lines.append(to_add)
    return lines

def load_patient_yaml(yaml_file):
    with open(yaml_file) as f:
        data = [x.split('\t') for x in yaml.load(f, Loader=yaml.SafeLoader)]
        return data

this_file_path = os.path.realpath(__file__)
sample_data_path = os.path.join(os.path.dirname(this_file_path), 'sample_data')

def test_run():
    hlas_file = os.path.join(sample_data_path, 'hla_data1.yaml')
    patients_file = os.path.join(sample_data_path, 'patient_data1.yaml')

    patients = PHAGE.getPatients(load_patient_yaml(patients_file), simple=1)
    hlas = load_hla_data(hlas_file)

    protein = 'Gag'

    epitopes_file = os.path.join(os.path.dirname(os.path.dirname(this_file_path)), 'epitopes.txt')

    epitopes = Epitope.Epitope.parseEpitopes(epitopes_file)
    print(epitopes[0])
    for e in epitopes:
        for i,hla in enumerate(e.hlas):
            e.hlas[i] = PHAGE.parseHLA(hla)
        for i,hla in enumerate(e.r2):
            e.r2[i] = PHAGE.parseHLA(hla)
        for i,hla in enumerate(e.r4):
            e.r4[i] = PHAGE.parseHLA(hla)
        e.hlas = set(e.hlas)
        e.r2 = set(e.r2)
        e.r4 = set(e.r4)

    print('hlas: {}'.format(hlas))
    for hla in hlas:
        hla[0] = PHAGE.parseHLA(hla[0])
    grouped_hlas = PHAGE.groupHLA(hlas)

    results = []
    done = set()
    for patient in patients:
        patient_results = PHAGE.analyzePatient(patients[patient], patient, protein, grouped_hlas, epitopes)
        results += patient_results

    results = sorted(results, key = lambda x: (x['pid'], x['pos'], x['hla']))

    PHAGE.displayResults(results, protein)