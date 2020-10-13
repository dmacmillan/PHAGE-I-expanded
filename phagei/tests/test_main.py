import os
import yaml
from phagei.Epitope import Epitope
import phagei.main as main

this_file_path = os.path.realpath(__file__)
sample_data_path = os.path.join(os.path.dirname(this_file_path), 'sample_data')


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


def test_run():
    hlas_file = os.path.join(sample_data_path, 'hla_data1.yaml')
    patients_file = os.path.join(sample_data_path, 'patient_data1.yaml')

    patients = main.getPatients(load_patient_yaml(patients_file), simple=1)
    hlas = load_hla_data(hlas_file)

    protein = 'Gag'

    epitopes_file = os.path.join(
        os.path.dirname(os.path.dirname(this_file_path)),
        'epitopes_v1.0.3.txt')

    epitopes = Epitope.parseEpitopes(epitopes_file)
    for e in epitopes:
        for i, hla in enumerate(e.hlas):
            e.hlas[i] = main.parseHLA(hla)
        for i, hla in enumerate(e.r2):
            e.r2[i] = main.parseHLA(hla)
        for i, hla in enumerate(e.r4):
            e.r4[i] = main.parseHLA(hla)
        e.hlas = set(e.hlas)
        e.r2 = set(e.r2)
        e.r4 = set(e.r4)

    for hla in hlas:
        hla[0] = main.parseHLA(hla[0])
    grouped_hlas = main.groupHLA(hlas)

    results = []
    done = set()
    for patient in patients:
        patient_results = main.analyzePatient(patients[patient], patient,
                                              protein, grouped_hlas, epitopes)
        results += patient_results

    results = sorted(results, key=lambda x: (x['pid'], x['pos'], x['hla']))
    # print(results)

    print(main.tsvResults(results, protein))


def test_gp41():
    hlas_file = os.path.join(sample_data_path, 'hla_data_gp41.yaml')
    patients_file = os.path.join(sample_data_path, 'patient_data_gp41.yaml')

    patients = main.getPatients(load_patient_yaml(patients_file), simple=1)
    hlas = load_hla_data(hlas_file)

    protein = 'gp41'

    epitopes_file = os.path.join(
        os.path.dirname(os.path.dirname(this_file_path)),
        'epitopes_v1.0.3.txt')

    epitopes = Epitope.parseEpitopes(epitopes_file)
    for e in epitopes:
        for i, hla in enumerate(e.hlas):
            e.hlas[i] = main.parseHLA(hla)
        for i, hla in enumerate(e.r2):
            e.r2[i] = main.parseHLA(hla)
        for i, hla in enumerate(e.r4):
            e.r4[i] = main.parseHLA(hla)
        e.hlas = set(e.hlas)
        e.r2 = set(e.r2)
        e.r4 = set(e.r4)

    for hla in hlas:
        hla[0] = main.parseHLA(hla[0])
    grouped_hlas = main.groupHLA(hlas)

    results = []
    done = set()
    for patient in patients:
        patient_results = main.analyzePatient(patients[patient], patient,
                                              protein, grouped_hlas, epitopes)
        results += patient_results

    print(main.tsvResults(results, protein))