# -*- coding: utf-8 -*-

import os
import sys

from datetime import datetime


def points_list_converted(points):
    """
    Converts list of curve points to match for
    a creating StudioSpline curve
    :param points: list
        List of curves points that has been obtained
        from dat file
    :return: list of coordinates corresponding to
        create_spline_from_points method of the NXAssembly class
    """

    coordinates = []
    for point in points:
        x = point[1]
        y = 0
        z = point[0]
        coordinates.append([x, y, z])
    return coordinates


class CommonClass:
    """
    A common class handles directories, files and etc.

    ...
    Methods
    -------
    platform_detect()
        Static method gets type of operating system
    new_file_name(file_name: str)
        Static method gets file name and returns
        new file name as string variable
    create_dir(dir_name: str)
        Static method creates dir
    get_files(path: str, ext: tuple)
        Static method returns files list with extension
        which is contained in tuple ext
    del_files(path: str, ext: tuple)
        Static method deletes files with extensions which are
        contained in tuple ext
    get_data_from_file(in_file=None, first_line=0)
        Static method reads file on the first_line and
        returns list of data

    """

    def __init__(self):
        super().__init__()

    @staticmethod
    def platform_detect():
        """
        :return: Current platform name, string type
        """
        if sys.platform.startswith('win'):
            current_platform = 'windows'
        elif sys.platform.startswith('linux'):
            current_platform = 'linux'
        elif sys.platform.startswith('aix'):
            current_platform = 'aix'
        elif sys.platform.startswith('cygwin'):
            current_platform = 'cigwin'
        elif sys.platform.startswith('darwin'):
            current_platform = 'macos'
        else:
            current_platform = 'unknown'
        return current_platform

    @staticmethod
    def new_file_name(file_name: str):
        """
        :param file_name: Old file name, type str, that will be renamed
        :return: str
            New file name, file name format old_file_name_HHMMSS_ddmmyy
        """
        path_, file_ = os.path.split(file_name)
        file_, ext_ = os.path.splitext(file_)

        dt = datetime.now()
        time_ = dt.strftime("%H%M%S")
        date_ = dt.strftime("%d%m%Y")
        new_file_name = file_ + "_" + date_ + "_" + time_ + ext_
        new_file_name = os.path.join(path_, new_file_name)
        return new_file_name

    @staticmethod
    def create_dir(dir_name: str):
        """
        :param dir_name: full directory name that will be created
        :return: Bool, str
            True if directory was created or directory exist, or False
            otherwise, and logging message
        """
        d = os.path.split(dir_name)
        try:
            os.mkdir(dir_name)
            msg = f"directory {d[1]} has been successfully created."
            return True, msg
        except PermissionError:
            msg = f"Directory '{d[1]}' has not been created. Access deny."
            return False, msg
        except FileExistsError:
            msg = f"Directory '{d[1]}' already exists."
            return True, msg
        finally:
            msg = f"Directory '{d[1]}' has not been created"
            return True, msg

    @staticmethod
    def get_files(path: str, ext: tuple):
        """
        :param path: Full directory name with extracting files
        :param ext: Tuple of extension of files that will be extracting
        :return: Bool, str
            True if files were got or False otherwise and logging message
        """

        arr_files = []
        try:
            with os.scandir(path) as dirs:
                for entry in dirs:
                    if os.path.splitext(entry.name)[1] in ext:
                        arr_files.append(os.path.join(path, entry.name))

            return True, arr_files
        except PermissionError:
            msg = "Iges files have not been obtained. "
            msg += f"Access denied to {os.path.split(path)[1]}."
            return False, msg
        except FileExistsError:
            msg = "prt files have not been created. "
            msg += f"Directory '{os.path.split(path)[1]}' doesn't exist."
            return False, msg
        except NotADirectoryError:
            msg = f"Trying to open directory '{path}' failed."
            msg += f"'{path}' is not a directory."
            return False, msg
        except FileNotFoundError:
            msg = f"Trying to open directory '{path}' failed. "
            msg += f"'{path}' doesn't exist."
            return False, msg

    @staticmethod
    def del_files(path: str, ext: tuple):
        """
        :param path: Full directory name which contains deleting files
        :param ext: Tuple of extensions of files which will be deleted
        :return: Bool, str
            True if files have been deleted or False otherwise and
                logging message
        """

        try:
            with os.scandir(path) as dirs:
                for entry in dirs:
                    if os.path.splitext(entry.name)[1] in ext:
                        os.remove(os.path.join(path, entry.name))
            msg = f"Redundant files have been deleted."
            return True, msg
        except FileExistsError:
            msg = "Removing excessive files. File doesn't exist."
            return False, msg
        except PermissionError:
            msg = "Removing excessive files. Access denied."
            return False, msg

    @staticmethod
    def get_data_from_file(in_file=None, first_line=0, sep=','):
        """
        :param in_file: Full file name from which data will be got
        :param first_line: First line from which data will be extracted
        :param sep: Separator for split method
        :return: Bool, str
            True if data was obtained or False otherwise and logging message

        """

        points = []
        if in_file:
            try:
                with open(in_file, 'r') as f:
                    content = f.readlines()
                    for point in content[first_line:]:
                        point = point.split(sep)
                        points.append(point)
                return True, points
            except PermissionError:
                msg = f"File {in_file} permission denied."
                return False, msg
            except FileExistsError:
                msg = f"File '{in_file}' doesn't exist."
                return False, msg
            except FileNotFoundError:
                msg = f"File '{in_file}' has not been found."
                return False, msg
