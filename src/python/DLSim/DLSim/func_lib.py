# -*- coding:utf-8 -*-
##############################################################
# Created Date: Friday, April 21st 2023
# Contact Info: luoxiangyong01@gmail.com
# Author/Copyright: Mr. Xiangyong Luo
##############################################################

import os
import datetime
from pathlib import Path
from typing import Union  # Python version <= 3.9
import ctypes
import platform


# A decorator to measure the time of a function
def func_running_time(func):
    def inner(*args, **kwargs):
        print(f'INFO Begin to run function: {func.__name__} …')
        time_start = datetime.datetime.now()
        res = func(*args, **kwargs)
        time_diff = datetime.datetime.now() - time_start
        print(
            f'INFO Finished running function: {func.__name__}, total: {time_diff.seconds}s')
        print()
        return res
    return inner


# convert OS path to standard linux path
def path2linux(path: Union[str, Path]) -> str:
    """Convert a path to a linux path, linux path can run in windows, linux and mac"""
    try:
        return path.replace("\\", "/")
    except Exception:
        return str(path).replace("\\", "/")


def get_filenames_from_folder_by_type(dir_name: str, file_type: str = "txt", isTraverseSubdirectory: bool = False) -> list:
    """Get all files in the folder with the specified file type

    Args:
        dir_name (str)                         : the folder path
        file_type (str, optional)              : the exact file type to specify, if file_type is "*" or "all", return all files in the folder. Defaults to "txt".
        isTraverseSubdirectory (bool, optional): get files inside the subfolder or not, if True, will traverse all subfolders. Defaults to False.

    Returns:
        list: a list of file paths
    """

    if isTraverseSubdirectory:
        files_list = []
        for root, dirs, files in os.walk(dir_name):
            files_list.extend([os.path.join(root, file) for file in files])
        if file_type in {"*", "all"}:
            return [path2linux(file) for file in files_list]
        return [path2linux(file) for file in files_list if file.split(".")[-1] == file_type]
    print("input dir:", dir_name, "\ninput file type", file_type)
    # files in the first layer of the folder
    if file_type in {"*", "all"}:
        return [path2linux(os.path.join(dir_name, file)) for file in os.listdir(dir_name)]
    return [path2linux(os.path.join(dir_name, file)) for file in os.listdir(dir_name) if file.split(".")[-1] == file_type]


def check_required_files_exist(required_files: list, dir_files: list) -> bool:
    # format the required file name to standard linux path
    required_files = [path2linux(os.path.abspath(filename)) for filename in required_files]

    required_files_short = [filename.split("/")[-1] for filename in required_files]
    dir_files_short = [filename.split("/")[-1] for filename in dir_files]

    # mask have the same length as required_files
    mask = [file in dir_files_short for file in required_files_short]
    if all(mask):
        return True

    print(f"Error: Required files are not satisfied, \
          missing files are: {[required_files_short[i] for i in range(len(required_files_short)) if not mask[i]]}")

    return False


def load_cpp_shared_library():
    _os = platform.system()
    if _os.startswith('Windows'):
        _dtalite_dll = os.path.join(os.path.dirname(__file__), 'pydtalite_bin/DTALite.dll')
    elif _os.startswith('Linux'):
        _dtalite_dll = os.path.join(os.path.dirname(__file__), 'pydtalite_bin/DTALite.so')
    elif _os.startswith('Darwin'):
        # check CPU is Intel or Apple Silicon
        if platform.machine().startswith('x86_64'):
            _dtalite_dll = os.path.join(os.path.dirname(__file__), 'pydtalite_bin/DTALite_x86.dylib')
        else:
            _dtalite_dll = os.path.join(os.path.dirname(__file__), 'pydtalite_bin/DTALite_arm.dylib')
    else:
        raise Exception('Please build the shared library compatible to your OS\
                        using source files')

    _dtalite_engine = ctypes.cdll.LoadLibrary(_dtalite_dll)

    _dtalite_engine.network_assignment.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int]
    return _dtalite_engine
