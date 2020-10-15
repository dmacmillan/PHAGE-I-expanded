import os
import yaml
import pytest
from tempfile import TemporaryFile
from pathlib import Path
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


@pytest.fixture
def setup():
    cwd = Path(os.path.realpath(__file__)).parent
    sample_data_path = cwd / 'sample_data'
    hlas_file = sample_data_path / 'hla_data1.yaml'
    patients_file = sample_data_path / 'patient_data1.yaml'

    hlas = load_hla_data(hlas_file)
    patients = main.getPatients(load_patient_yaml(patients_file), simple=1)
    protein = 'Gag'

    hlas_tempfile = './tmp_hlas'
    with open(hlas_tempfile, 'w') as o:
        o.write('\n'.join(['\t'.join(x) for x in hlas]))
    patients_tempfile = './tmp_patients'
    with open(patients_tempfile, 'w') as o:
        o.write('\n'.join(['\t'.join(x) for x in patients]))

    return {'cwd': cwd, 'patients': patients_tempfile, 'hlas': hlas_tempfile}


def test1(setup):
    tsv_results = main.run(setup['hlas'], setup['patients'], 'Gag')
    print(tsv_results)
    os.remove(setup['hlas'])
    os.remove(setup['patients'])


# def test_run():
#     hlas_file = os.path.join(sample_data_path, 'hla_data1.yaml')
#     patients_file = os.path.join(sample_data_path, 'patient_data1.yaml')

#     patients = main.getPatients(load_patient_yaml(patients_file), simple=1)
#     hlas = load_hla_data(hlas_file)

#     protein = 'Gag'

#     epitopes_file = os.path.join(
#         os.path.dirname(os.path.dirname(this_file_path)),
#         'epitopes_v1.0.3.txt')

#     epitopes = Epitope.parseEpitopes(epitopes_file)
#     for e in epitopes:
#         for i, hla in enumerate(e.hlas):
#             e.hlas[i] = main.parseHLA(hla)
#         for i, hla in enumerate(e.r2):
#             e.r2[i] = main.parseHLA(hla)
#         for i, hla in enumerate(e.r4):
#             e.r4[i] = main.parseHLA(hla)
#         e.hlas = set(e.hlas)
#         e.r2 = set(e.r2)
#         e.r4 = set(e.r4)

#     for hla in hlas:
#         hla[0] = main.parseHLA(hla[0])
#     grouped_hlas = main.groupHLA(hlas)

#     results = []
#     done = set()
#     for patient in patients:
#         patient_results = main.analyzePatient(patients[patient], patient,
#                                               protein, grouped_hlas, epitopes)
#         results += patient_results

#     results = sorted(results, key=lambda x: (x['pid'], x['pos'], x['hla']))
#     # print(results)

#     print(main.tsvResults(results, protein))

# def test_gp41():
#     hlas_file = os.path.join(sample_data_path, 'hla_data_gp41.yaml')
#     patients_file = os.path.join(sample_data_path, 'patient_data_gp41.yaml')

#     patients = main.getPatients(load_patient_yaml(patients_file), simple=1)
#     hlas = load_hla_data(hlas_file)

#     protein = 'gp41'

#     epitopes_file = os.path.join(
#         os.path.dirname(os.path.dirname(this_file_path)),
#         'epitopes_v1.0.3.txt')

#     epitopes = Epitope.parseEpitopes(epitopes_file)
#     for e in epitopes:
#         for i, hla in enumerate(e.hlas):
#             e.hlas[i] = main.parseHLA(hla)
#         for i, hla in enumerate(e.r2):
#             e.r2[i] = main.parseHLA(hla)
#         for i, hla in enumerate(e.r4):
#             e.r4[i] = main.parseHLA(hla)
#         e.hlas = set(e.hlas)
#         e.r2 = set(e.r2)
#         e.r4 = set(e.r4)

#     for hla in hlas:
#         hla[0] = main.parseHLA(hla[0])
#     grouped_hlas = main.groupHLA(hlas)

#     results = []
#     done = set()
#     for patient in patients:
#         patient_results = main.analyzePatient(patients[patient], patient,
#                                               protein, grouped_hlas, epitopes)
#         results += patient_results

#     print(main.tsvResults(results, protein))