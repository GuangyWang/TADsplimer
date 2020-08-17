import numpy as np
import os
import cv2
import math


def TAD_calling(file, up, down):
    #    warnings.filterwarnings("ignore", category=RRuntimeWarning)
    if os.stat(file).st_size != 0:
        os.system(
            'Rscript '+os.path.join(os.path.dirname(os.path.realpath(__file__)), 'TAD_calling.R') +
            ' -f ' + file + ' -u ' + str(up) + ' -d ' + str(down) +
            ' -o ' + file + '.boundary.txt')
        TAD = np.loadtxt(file + '.boundary.txt')
        TAD = TAD[:-1, :]
    else:
        TAD = np.array([0, 0, 0])
    return(TAD)


def file_split(file, outfile, length=1000, print_subcontact=0, aliases='None'):
    '''
    divide whole chromosome into different sub chromosome files and calculate cutoffs
    :param file: input chromosome file
    :param outfile: output path for sub chromosome files
    :param length: length of sub chromosome
    :param print_subcontact: if print_subcontact==1, not save sub chromosome files, if print_subcontact!=1, save sub
    chromosome files
    :return: up limit and down limit for this chromosome
    '''
    base = os.path.basename(file)
    if aliases != "None":
        filename = aliases
    else:
        filename = os.path.splitext(base)[0]

    print(file)
    contact_map = np.loadtxt(file)
    r, c = contact_map.shape
    num = int(r/length)
    mean_list = []
    files = []

    if num>0:
        flag = 0
        if (num*length+length/2) > c:
            flag = 1

        if flag == 1:
            for i in range(num):
                sub_contact = contact_map[(i * length):((i + 1) * length), (i * length):((i + 1) * length)]
                out = os.path.join(outfile, filename + '.' + str(i * length) + '.' + str((i + 1) * length) + '.subchr')
                files.append(out)
                mean_list.append(np.mean(sub_contact))
                if print_subcontact != 1:
                    np.savetxt(out, sub_contact, delimiter='\t')
                if i != (num-1):
                    sub_contact2 = contact_map[int(i * length+length/2):int((i + 1) * length+length/2),
                                   int(i * length+length/2):int((i + 1) * length+length/2)]
                    out2 = os.path.join(outfile,
                                       filename + '.' + str(int(i * length+length/2)) + '.' + str(int((i + 1) * length + length/2)) + '.subchr')
                    files.append(out2)
                    mean_list.append(np.mean(sub_contact2))
                    if print_subcontact != 1:
                        np.savetxt(out2, sub_contact2, delimiter='\t')
            sub_contact = contact_map[(num * length):c, (num * length):c]
            out = os.path.join(outfile, filename + '.' + str(num * length) + '.' + str(c) + '.subchr')
            files.append(out)
            mean_list.append(np.mean(sub_contact))
            if print_subcontact != 1:
                np.savetxt(out, sub_contact, delimiter='\t')
            sub_contact2 = contact_map[int((num-1) * length + length / 2):int(((num-1) + 1) * length + length / 2),
                           int((num-1) * length + length / 2):int(((num-1) + 1) * length + length / 2)]
            out2 = os.path.join(outfile, filename + '.' + str(int((num-1) * length + length / 2)) + '.' + str(int(((num-1) + 1) * length + length / 2)) + '.subchr')
            files.append(out2)
            mean_list.append(np.mean(sub_contact2))
            if print_subcontact != 1:
                np.savetxt(out2, sub_contact2, delimiter='\t')
        else:
            for i in range(num):
                sub_contact = contact_map[(i * length):((i + 1) * length), (i * length):((i + 1) * length)]
                out = os.path.join(outfile, filename + '.' + str(i * length) + '.' + str((i + 1) * length) + '.subchr')
                files.append(out)
                mean_list.append(np.mean(sub_contact))
                if print_subcontact != 1:
                    np.savetxt(out, sub_contact, delimiter='\t')
                sub_contact2 = contact_map[int(i * length+length/2):int((i + 1) * length+length/2),
                               int(i * length+length/2):int((i + 1) * length+length/2)]
                out2 = os.path.join(outfile,
                                   filename + '.' + str(int(i * length+length/2)) + '.' + str(int((i + 1) * length + length/2)) + '.subchr')
                files.append(out2)
                mean_list.append(np.mean(sub_contact2))
                if print_subcontact != 1:
                    np.savetxt(out2, sub_contact2, delimiter='\t')
            sub_contact = contact_map[(num * length):c, (num * length):c]
            out = os.path.join(outfile, filename + '.' + str(num * length) + '.' + str(c) + '.subchr')
            files.append(out)
            mean_list.append(np.mean(sub_contact))
            if print_subcontact != 1:
                np.savetxt(out, sub_contact, delimiter='\t')
            sub_contact2 = contact_map[int(num * length + length / 2):int((num + 1) * length + length / 2),
                           int(num * length + length / 2):int((num + 1) * length + length / 2)]
            out2 = os.path.join(outfile,
                                filename + '.' + str(int(num * length + length / 2)) + '.' + str(int(
                                    (num + 1) * length + length / 2)) + '.subchr')
            files.append(out2)
            mean_list.append(np.mean(sub_contact2))
            if print_subcontact != 1:
                np.savetxt(out2, sub_contact2, delimiter='\t')
    else:
        out = os.path.join(outfile, filename + '.' + str(0) + '.' + str(r) + '.subchr')
        files.append(out)
        if print_subcontact != 1:
            np.savetxt(out, contact_map, delimiter='\t')
    mean_list = np.array(mean_list)
    mean_list = mean_list[~np.isnan(mean_list)]
    cutoff_up = np.median(np.array(mean_list) / 0.22)
    cutoff_down = np.median(np.array(mean_list) / 2)
    files.sort()
    return cutoff_up, cutoff_down, files


def TAD_plot(f1, output, TAD, down, up):
    '''
    :param f1: contact map
    :param down: down limit for color key
    :param up: up limit for color key
    :return:
    '''
    if os.stat(f1).st_size != 0:
        img = np.loadtxt(f1)
        img[img > up] = up
        img[img < down] = 0
        img = img * (255 / up)
        img = img.astype(np.uint8)
        color_img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
        color_img[:, :, 0] = 255 - color_img[:, :, 0]
        color_img[:, :, 1] = 255 - color_img[:, :, 1]
        color_img[:, :, 2] = 255
        for x in TAD:
            if float(x[0]) > 0.8:
                color_img = cv2.rectangle(color_img, (math.ceil(float(x[1])), math.ceil(float(x[1]))),
                                          (math.ceil(float(x[2])), math.ceil(float(x[2]))), (0, 0, 0), 2,
                                          lineType=4)
                cv2.imwrite(output, color_img)
    else:
        print('file '+f1+' is empty')



def subTAD_calling(path, up, down, plt=1, prnt=1, aliase='None'):
    dirs = os.listdir(path)
    tad = np.empty(shape=[0, 3])
    for file in dirs:
        if file[-7:] == '.subchr':

            print(file)

            start = file.split('.')[-3]
            tad_temp = TAD_calling(os.path.join(path, file), up, down)
            tad_location = tad_temp

            if plt == 1:
                try:
                    TAD_plot(os.path.join(path, file), os.path.join(path, file + '.tiff'),
                                         tad_temp, down, up)
                except:
                    pass
            try:
                tad_location[:, 1:3] = tad_location[:, 1:3] + int(start)
                if not np.array_equal(tad_location, np.array([0, int(start), int(start)])):
                    tad = np.vstack((tad, tad_location))
                if prnt == 1:
                    np.savetxt(os.path.join(path, file + '.TAD.txt'), tad_location, delimiter="\t")
            except:
                tad_location = np.array([[0, 0, 0]])
                if prnt == 1:
                    np.savetxt(os.path.join(path, file + '.TAD.txt'), tad_location, delimiter="\t")
    sorted_idx = np.lexsort(tad[:, 1:3].T)
    sorted_tad = tad[sorted_idx, :]
    row_mask = np.append([True], np.any(np.diff(sorted_tad[:, 1:3], axis=0), 1))
    if sorted_tad.size != 0:
        unique_tad = sorted_tad[row_mask]
        if aliase == 'None':
            np.savetxt(os.path.join(path, 'TAD.merge.txt'), unique_tad, delimiter="\t")
        else:
            np.savetxt(os.path.join(path, aliase + '.TAD.merge.txt'), unique_tad, delimiter="\t")
    else:
        print('No split/merged TAD\n')
