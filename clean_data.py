import numpy as np


def txt_to_csv(data_file, out_file):
    results_arr = np.loadtxt(data_file)
    print(results_arr.shape)
    x_data = results_arr[0, :]
    mean_data = np.mean(results_arr[1:, :], axis=0)
    std_data = np.std(results_arr[1:, :], axis=0)
    print(len(x_data), len(mean_data), len(std_data))
    with open(out_file, 'w+') as f:
        for i in range(len(x_data)):
            f.write(str(x_data[i]) + ',' + str(mean_data[i]) + ',' + str(std_data[i]) + '\n')


def txt_to_csv_1d(data_file, out_file):
    data = np.loadtxt(data_file)
    print(data.shape)
    col1 = data[0,:]
    col2 = data[1,:]
    assert len(col1) == len(col2), 'columns are not of a same size!'
    with open(out_file, 'w+') as f:
        for i in range(len(col1)):
            f.write(str(col1[i]) + ',' + str(col2[i]) + '\n')



data_file = 'surface_roughness\\sf_20_50_20_current_rough.txt'
out_file = 'surface_roughness\\sf_20_50_20_current_ideal.csv'
txt_to_csv_1d(data_file, out_file)
