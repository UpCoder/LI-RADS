# -*- coding: utf-8 -*-
import numpy as np
import pydicom
import os
import SimpleITK as itk
from PIL import Image, ImageDraw
import Queue
import gc
import copy


# 调整窗宽 窗位
def rejust_pixel_value(image):
    image = np.array(image)
    ww = np.float64(250)
    wc = np.float64(55)
    ww = max(1, ww)
    lut_min = 0
    lut_max = 255
    lut_range = np.float64(lut_max) - lut_min

    minval = wc - ww / 2.0
    maxval = wc + ww / 2.0
    image[image < minval] = minval
    image[image > maxval] = maxval
    to_scale = (minval <= image) & (image <= maxval)
    image[to_scale] = ((image[to_scale] - minval) / (ww * 1.0)) * lut_range + lut_min
    return image


# 读取单个DICOM文件
def read_file(file_name):
    header = pydicom.read_file(file_name)
    image = header.pixel_array
    image = header.RescaleSlope * image + header.RescaleIntercept
    # print header.RescaleSlope,header.RescaleIntercept

    return image


# 读取DICOM文件序列
def read_dicom_series(dir_name):
    print 'read dicom ', dir_name
    files = list(os.listdir(dir_name))
    files.sort()
    res = []
    for file in files:
        if file.endswith('DCM'):
            cur_file = os.path.join(dir_name, file)
            res.append(read_file(cur_file))
    res = np.array(res)
    return res

def read_dicom_series_itk(dir_name):
    reader=itk.ImageSeriesReader()
    dicom_series=reader.GetGDCMSeriesFileNames(dir_name)
    reader.SetFileNames(dicom_series)
    images=reader.Execute()
    image_array = itk.GetArrayFromImage(images)
    return image_array

def dicom_series_mhd(dir_name, save_path):
    reader = itk.ImageSeriesReader()
    dicom_series = reader.GetGDCMSeriesFileNames(dir_name)
    reader.SetFileNames(dicom_series)
    images = reader.Execute()
    itk.WriteImage(images, save_path)



# 读取mhd文件
def read_mhd_image(file_path):
    header = itk.ReadImage(file_path)
    image = itk.GetArrayFromImage(header)
    return np.array(image)


# 保存mhd文件
def save_mhd_image(image, file_name):
    print 'image type is ', type(image)
    header = itk.GetImageFromArray(image)
    itk.WriteImage(header, file_name)


# 一次读取多个ｍｈｄ文件，然后统一放缩至指定大小
def read_mhd_images(paths, new_size=None, avg_liver_values=None):
    images = []
    for index, path in enumerate(paths):
        # print path
        cur_image = read_mhd_image(path)
        cur_img = np.asarray(cur_image, np.float32)
        if avg_liver_values is not None:
            for i in range(0):
                cur_img = cur_img * cur_img
                cur_img = cur_img / avg_liver_values[index]
        if new_size is not None:
            cur_img = Image.fromarray(np.asarray(cur_img, np.float32))
            print new_size
            cur_img = cur_img.resize(new_size)
            #print np.shape(cur_img), path
            cur_image = np.array(cur_img)
        images.append(cur_image)
    return images


# 将灰度图像转化为RGB通道
def conver_image_RGB(gray_image):
    shape = list(np.shape(gray_image))
    image_arr_rgb = np.zeros(shape=[shape[0], shape[1], 3])
    image_arr_rgb[:, :, 0] = gray_image
    image_arr_rgb[:, :, 1] = gray_image
    image_arr_rgb[:, :, 2] = gray_image
    return image_arr_rgb


# 将一个矩阵保存为图片
def save_image(image_arr, save_path):
    if len(np.shape(image_arr)) == 2:
        image_arr = conver_image_RGB(image_arr)
    image = Image.fromarray(np.asarray(image_arr, np.uint8))
    image.save(save_path)


# 将图像画出来，并且画出标记的病灶
def save_image_with_mask(image_arr, mask_image, save_path):
    shape = list(np.shape(image_arr))
    image_arr_rgb = np.zeros(shape=[shape[0], shape[1], 3])
    image_arr_rgb[:, :, 0] = image_arr
    image_arr_rgb[:, :, 1] = image_arr
    image_arr_rgb[:, :, 2] = image_arr
    image = Image.fromarray(np.asarray(image_arr_rgb, np.uint8))
    image_draw = ImageDraw.Draw(image)
    [ys, xs] = np.where(mask_image != 0)
    miny = np.min(ys)
    maxy = np.max(ys)
    minx = np.min(xs)
    maxx = np.max(xs)
    ROI = image_arr_rgb[miny-1:maxy+1, minx-1:maxx+1, :]
    ROI_Image = Image.fromarray(np.asarray(ROI, np.uint8))

    for index, y in enumerate(ys):
        image_draw.point([xs[index], y], fill=(255, 0, 0))
    image.save(save_path)
    ROI_Image.save(os.path.join(os.path.dirname(save_path), os.path.basename(save_path).split('.')[0]+'_ROI.jpg'))
    del image, ROI_Image
    gc.collect()
# 获取单位方向的坐标，
# 比如dim=2，则返回的数组就是[[-1, -1],...[1, 1]]
def get_direction_index(dim=3, cur_dir=[]):
    res = []
    for i in range(-1, 2):
        cur_dir.append(i)
        if dim != 1:
            res.extend(get_direction_index(dim-1, cur_dir))
            cur_dir.pop()
        else:
            res.append(copy.copy(cur_dir))
            cur_dir.pop()
    return res


# 验证数组arr1的值是否在arr_top 和 arr_dowm之间
# arr_top[i] > arr1[i] >= arr_down[i]
def value_valid(arr1, arr_top, arr_down):
    for index, item in enumerate(arr1):
        if arr_down[index] <= item < arr_top[index]:
            continue
        else:
            return False
    return True


# 将mask文件中的多个病灶拆分出来
def split_mask_image(total_mask_image_path, save_paths):
    directions = get_direction_index(dim=3)

    def find_connected_components(position, mask_image, flag):
        queue = Queue.Queue()
        points = []
        queue.put(position)
        while not queue.empty():
            cur_position = queue.get()
            points.append(cur_position)
            for direction in directions:
                new_z = cur_position[0] + direction[0]
                new_y = cur_position[1] + direction[1]
                new_x = cur_position[2] + direction[2]
                if value_valid([new_z, new_y, new_x], np.shape(mask_image), [0, 0, 0])\
                        and flag[new_z, new_y, new_x] == 1 \
                        and mask_image[new_z, new_y, new_x] != 0:
                    queue.put([new_z, new_y, new_x])
                    flag[new_z, new_y, new_x] = 0
        return points
    mask_image = read_mhd_image(total_mask_image_path)
    mask_image = np.array(mask_image)
    [z, y, x] = np.shape(mask_image)
    flag = np.ones(
        shape=[
            z, y, x
        ]
    )
    flag[np.where(mask_image == 0)] = 0
    index = 0
    for o in range(z):
        for n in range(y):
            for m in range(x):
                if mask_image[o, n, m] != 0 and flag[o, n, m] == 1:
                    flag[o, n, m] = 0
                    points = find_connected_components([o, n, m], mask_image, flag)
                    new_mask = np.zeros(
                        shape=[z, y, x]
                    )
                    for point in points:
                        new_mask[point[0], point[1], point[2]] = 1
                    save_mhd_image(new_mask, save_paths[index])
                    index += 1
                    print len(points)






# 获取label的分布
def get_distribution_label(labels):
    min_value = np.min(labels)
    max_value = np.max(labels)
    my_dict = {}
    for label in labels:
        if label in my_dict:
            my_dict[label] += 1
        else:
            my_dict[label] = 0
    return my_dict


#　将数据打乱
def shuffle_image_label(images, labels):
    labels = np.array(labels)
    if len(images) == 2:
        random_index = range(len(images[0]))
    else:
        random_index = range(len(images))
    np.random.shuffle(random_index)
    labels = labels[random_index]
    new_images = []
    for cur_index in random_index:
        new_images.append(images[cur_index])
    return new_images, labels


# 将数据按照指定的方式排序
def changed_shape(image, shape):
    new_image = np.zeros(
        shape=shape
    )
    batch_size = shape[0]
    for z in range(batch_size):
        for phase in range(shape[-1]):
            if shape[-1] == 1:
                new_image[z, :, :, phase] = image[z]
            else:
                new_image[z, :, :, phase] = image[z, phase]
    del image
    gc.collect()
    return new_image


# 将一幅图像的mask外部全部标记为０
def mark_outer_zero(image, mask_image):
    def is_in(x, y, mask):
        sum1 = np.sum(mask[0:x, y])
        sum2 = np.sum(mask[x:, y])
        sum3 = np.sum(mask[x, 0:y])
        sum4 = np.sum(mask[x, y:])
        if sum1 != 0 and sum2 != 0 and sum3 != 0 and sum4 != 0:
            return True
        return False
    def fill_region(mask, x, y):
        queue = Queue.Queue()
        queue.put([x, y])
        directions = [[1, 0], [-1, 0], [0, 1], [0, -1]]
        while not queue.empty():
            point = queue.get()
            for direction in directions:
                new_x = point[0] + direction[0]
                new_y = point[1] + direction[1]
                if value_valid([new_x, new_y], np.shape(mask), [0, 0]) and mask[new_x, new_y] == 0:
                    mask[new_x, new_y] = 1
                    queue.put([new_x, new_y])
        return mask
    [w, h] = np.shape(image)
    mask_image_copy = mask_image.copy()
    for i in range(w):
        # find = False
        for j in range(h):
            if mask_image_copy[i, j] == 0:
                if is_in(i, j, mask_image_copy):
                    print i, j
                    mask_image_copy[i, j] = 1
                    mask_image_copy = fill_region(mask_image_copy, i, j)
                    # find = True
                    image[np.where(mask_image_copy == 0)] = 0
                    return image, mask_image_copy
    print 'Error'
    return image, mask_image_copy


# 显示一幅图像
def show_image(image_arr, title=None):
    image = Image.fromarray(image_arr)
    image.show(title=title)


# 计算针对二分类错了多少个
def acc_binary_acc(logits, label):
    acc_count = 0.0
    logits = copy.copy(logits)
    label = copy.copy(label)
    logits = np.array(logits)
    label = np.array(label)
    logits[logits == 1] = 0
    logits[logits == 3] = 0
    logits[logits == 2] = 1
    logits[logits == 4] = 1

    label[label == 1] = 0
    label[label == 3] = 0
    label[label == 2] = 1
    label[label == 4] = 1
    for index, logit in enumerate(logits):
        if label[index] == logit:
            acc_count += 1
    return (1.0 * acc_count) / (1.0 * len(logits))

# 计算Ａｃｃｕｒａｃｙ，并且返回每一类最大错了多少个
def calculate_acc_error(logits, label, show=True):
    error_index = []
    error_dict = {}
    error_dict_record = {}
    error_count = 0
    error_record = []
    label = np.array(label)
    for index, logit in enumerate(logits):
        if logit != label[index]:
            error_count += 1
            if label[index] in error_dict.keys():
                error_dict[label[index]] += 1   # 该类别分类错误的个数加１
                error_dict_record[label[index]].append(logit)   # 记录错误的结果
            else:
                error_dict[label[index]] = 1
                error_dict_record[label[index]] = [logit]
            error_index.append(index)
            error_record.append(logit)
    acc = (1.0 * error_count) / (1.0 * len(label))
    if show:
        for key in error_dict.keys():
            print 'label is %d, error number is %d, all number is %d, acc is %g'\
                  % (key, error_dict[key], np.sum(label == key), 1-(error_dict[key]*1.0)/(np.sum(label == key) * 1.0))
            print 'error record　is ', error_dict_record[key]
    return error_dict, error_dict_record, acc, error_index, error_record


def get_shuffle_index(n):
    random_index = range(n)
    np.random.shuffle(random_index)
    return random_index


# 对所有病灶进行线性增强
def linear_enhancement(path='/home/give/PycharmProjects/MedicalImage/imgs/LiverAndLesion_bg/TRAIN'):
    image_dirs = os.listdir(path)
    for image_dir in image_dirs:
        for phase in ['nc', 'art', 'pv']:
            liver_path = os.path.join(path, image_dir, 'liver_' + phase + '.jpg')
            tumor_path = os.path.join(path, image_dir, 'tumor_' + phase + '.jpg')
            liver_image = Image.open(liver_path)
            tumor_image = Image.open(tumor_path)
            from Slice.TwoFolderMaxSlice.Slice_Base_Liver_Tumor import Liver_Tumor_Operations
            new_tumor = Liver_Tumor_Operations.tumor_linear_enhancement_only_tumor(tumor_image=np.array(tumor_image))
            new_path = os.path.join(path, image_dir, 'tumor_' + phase + '_enhancement.jpg')
            image = Image.fromarray(np.asarray(new_tumor, np.uint8))
            image.save(new_path)


def extract_avg_liver_dict(txt_path='/home/give/Documents/dataset/MedicalImage/MedicalImage/average_pixel_value.txt'):
    filed = open(txt_path)
    lines = filed.readlines()
    res_dict = {}
    for line in lines:
        split_res = line.split(',')
        srrid = split_res[0]
        srrid = int(srrid.split('-')[0])
        avg_liver = [float(split_res[i]) for i in range(1, 4)]
        res_dict[srrid] = avg_liver
    return res_dict


def calculate_tp(logits, labels):
    count = 0
    for index, logit in enumerate(logits):
        if logit == labels[index] and logit == 1:
            count += 1
    return count

def calculate_recall(logits, labels):
    tp = calculate_tp(logits=logits, labels=labels)
    recall = (tp * 1.0) / (np.sum(labels == 1) * 1.0)
    return recall


def calculate_precision(logits, labels):
    tp = calculate_tp(logits=logits, labels=labels)
    precision = (tp * 1.0) / (np.sum(logits == 1) * 1.0)
    return precision


def get_game_evaluate(logits, labels, argmax=None):
    logits = np.array(logits)
    labels = np.array(labels)
    if argmax is not None:
        logits = np.argmax(logits, argmax)
        labels = np.argmax(labels, argmax)
    recall = calculate_recall(logits=logits, labels=labels)
    precision = calculate_precision(logits=logits, labels=labels)
    f1_score = (2*precision*recall) / (precision + recall)
    return recall, precision, f1_score

if __name__ == '__main__':
    linear_enhancement()