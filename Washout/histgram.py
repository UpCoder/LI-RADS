# -*- coding=utf-8 -*-
from Tools.Tools import *
from matplotlib import pyplot as plt
from glob import glob
from Config import WashOutConfig
import xlrd
from skimage.morphology import erosion, square


def get_histgram(mhd_image, mask_image, liver_mask, save_path=None):
    '''
    根据ｍｈｄ文件和ｍａｓｋ文件得到病灶的像素值分布和肝脏的像素值分布
    :param mhd_image:
    :param mask_image: mask value = 1 represent the pixel is lesion
    :param liver_mask: mask value =1 represent the pixel is liver
    :return: 255 * 1 的数组
    '''
    def hexencode(rgb):
        r = rgb[0]
        g = rgb[1]
        b = rgb[2]
        a = rgb[3]
        return '#%02x%02x%02x%02x' % (r, g, b, a)
    xs ,ys = np.where(mask_image == 1)
    for index, x in enumerate(xs):
        y = ys[index]
        print x, y, mhd_image[x, y]
    lesion_values = mhd_image[mask_image == 1]
    liver_values = mhd_image[liver_mask == 1]
    lesion_values = lesion_values.flatten()
    liver_values = liver_values.flatten()
    liver_mean = np.mean(liver_values)

    for i in range(0):
        lesion_values = lesion_values * lesion_values / liver_mean
    print 'lesion max pixel is %d, liver max pixel is %d' % (np.max(lesion_values), np.max(liver_values))
    print 'lesion mean pixel is %d, liver mean pixel is %d' % (np.mean(lesion_values), np.mean(liver_values))
    np.random.shuffle(lesion_values)
    np.random.shuffle(liver_values)
    min_pixel_num = min(len(lesion_values), len(liver_values))
    lesion_values = lesion_values[:min_pixel_num]
    liver_values = liver_values[:min_pixel_num]
    lesion_values[lesion_values < 0] = 0
    liver_values[liver_values < 0] = 0
    pixel_values = np.concatenate([lesion_values, liver_values], axis=0)
    (n, bins) = np.histogram(pixel_values, bins=255)
    (n_lesion, _) = np.histogram(lesion_values, bins=255)
    (n_liver, _) = np.histogram(liver_values, bins=255)
    print np.sum(n), min_pixel_num * 2
    n[n <5 ] = 0    # 去除离散点的影响
    n_liver[n_liver < 5] = 0
    n_lesion[n_lesion < 5] = 0
    # print zip(range(255), n_lesion)
    # plt.bar(range(255), n_liver, 1, color=hexencode([0, 255, 0, 128]), linewidth=2)
    # plt.bar(range(255), n_lesion, 1, color=hexencode([255, 0, 0, 128]), linewidth=2)
    # plt.xlabel("X-axis")
    # plt.ylabel("Y-axis")
    # plt.title("histgram")
    # plt.show()

    plt.hist(lesion_values, bins=255, color=hexencode([255, 0, 0, 128]))
    plt.hist(liver_values, bins=255, color=hexencode([0, 255, 0, 128]))

    if save_path is not None:
        plt.savefig(save_path)
    else:
        plt.show()
    plt.close()
    return n, n_lesion, n_liver

'''
    针对RSNA的数据生成直方图
'''
def generate_histgrams(RSNA_DIR, excel_path, save_dirs, rewrite=True):
    def process_float(float_val):
        str_val = str(float_val)
        str_val = str_val[:str_val.find('.')]
        return str_val
    tables = xlrd.open_workbook(excel_path)
    table = tables.sheets()[0]
    nRows = table.nrows
    label_dict = {}
    for r_index in range(1, nRows):
        col_values = table.row_values(r_index)
        pcid = process_float(col_values[2]) + '-' + process_float(col_values[3])
        if pcid in label_dict.keys():
            print 'Error pcid repeat\n'
            return
        label_dict[pcid] = int(col_values[10])
    sub_dirs = ['HCC_first', 'HCC_second']
    # sub_dirs = ['HCC_second']
    features = []
    splited_features = []
    labels = []
    for sub_dir in sub_dirs:
        cur_dir = os.path.join(RSNA_DIR, sub_dir)
        for dir_name in os.listdir(cur_dir):
            # 如果是relationtxt或者是删除的文件的话，跳过
            if not os.path.isdir(os.path.join(cur_dir, dir_name)) or dir_name in ['3892299-3594663']:
                continue
            cur_label = label_dict[dir_name]
            if not rewrite and os.path.exists(os.path.join(save_dirs[cur_label], dir_name+'.png')):
                continue
            pv_dicom_dir = os.path.join(cur_dir, dir_name, 'PV')
            pv_mask_path = glob(os.path.join(cur_dir, dir_name, '*_PV.mhd'))[0]
            pv_images = read_dicom_series(pv_dicom_dir)
            pv_masks = read_mhd_image(pv_mask_path)
            [zs, _, _] = np.where(pv_masks != 0)
            if dir_name in ['3715982-3405713']:
                mhd_image = pv_images[zs[0], :, :]
            else:
                mhd_image = pv_images[len(pv_images) - zs[0], :, :]
            mask_image = pv_masks[zs[0], :, :]
            liver_mask = copy.copy(mask_image)
            mask_image[mask_image != 1] = 0
            liver_mask[liver_mask == 1] = 0
            liver_mask[liver_mask == 2] = 1
            kernel_size = 11
            if np.sum(erosion(mask_image, square(kernel_size))) != 0:
                mask_image = erosion(mask_image, square(kernel_size))
                liver_mask = erosion(liver_mask, square(kernel_size))
            feature, feature_lesion, feature_liver = get_histgram(mhd_image, mask_image, liver_mask, os.path.join(save_dirs[cur_label], dir_name+'.png'))
            features.append(feature)
            splited_features.append(np.concatenate([feature_lesion, feature_liver], axis=0))
            labels.append(cur_label)
    np.save('./histgram_splited_features', splited_features)
    np.save('./histgram_features', features)
    np.save('./histgram_labels', labels)

def get_histgram_copy(mhd_image, mask_image, liver_mask, save_path=None):
    '''
    根据ｍｈｄ文件和ｍａｓｋ文件得到病灶的像素值分布和肝脏的像素值分布
    :param mhd_image:
    :param mask_image: mask value = 1 represent the pixel is lesion
    :param liver_mask: mask value =1 represent the pixel is liver
    :return: 255 * 1 的数组
    '''

    def hexencode(rgb):
        r = rgb[0]
        g = rgb[1]
        b = rgb[2]
        a = rgb[3]
        return '#%02x%02x%02x%02x' % (r, g, b, a)

    xs, ys = np.where(mask_image == 1)
    for index, x in enumerate(xs):
        y = ys[index]
        print x, y, mhd_image[x, y]
    # lesion_values = mhd_image[mask_image == 1]
    print mhd_image[0, 0]

if __name__ == '__main__':
    def show_care(image, mask_image):
        image[mask_image == 0] = 0
        return image
    # generate_histgrams(
    #     WashOutConfig.RSNA_DATA_DIR,
    #     '/home/give/Documents/dataset/LI-RADS/data/RSNA/RSNA.xlsx',
    #     [
    #         '/home/give/Documents/dataset/LI-RADS/data/RSNA/histgrams/washout/negative',
    #         '/home/give/Documents/dataset/LI-RADS/data/RSNA/histgrams/washout/positive'
    #     ]
    # )
    pcid = '3631930-3242626'
    id = '43'
    # pcid = '3950729-3693551'
    # id = '48'
    pv_dicom_dir = '/home/give/Documents/dataset/LI-RADS/data/RSNA/HCC_first/' + pcid + '/PV'
    pv_mask_path = '/home/give/Documents/dataset/LI-RADS/data/RSNA/HCC_first/' + pcid +'/' + id +'_PV.mhd'
    pv_images = read_dicom_series_itk(pv_dicom_dir)
    dicom_series_mhd(pv_dicom_dir, './pv.mhd')
    print np.shape(pv_images)

    pv_masks = read_mhd_image(pv_mask_path)
    [zs, ys, xs] = np.where(pv_masks != 0)
    print zs
    mhd_image = pv_images[len(pv_images) - zs[0], :, :]
    mask_image = pv_masks[zs[0], :, :]
    liver_mask = copy.copy(mask_image)
    mask_image[mask_image != 1] = 0
    liver_mask[liver_mask == 1] = 0
    liver_mask[liver_mask == 2] = 1
    from skimage.morphology import erosion, square
    # mask_image = erosion(mask_image, square(11))
    # liver_mask = erosion(liver_mask, square(11))


    get_histgram_copy(mhd_image, mask_image, liver_mask, None)