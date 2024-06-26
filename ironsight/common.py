# -*- coding: utf-8 -*-

# ~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import sys
import os
import inspect
from collections import *
from glob import iglob
import datetime
import logging
import json
import gzip
import subprocess
import stat

# Local imports
from ironsight import __version__ as pkg_version
from ironsight import __name__ as pkg_name

# Third party imports
import colorlog

# ~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~#
def opt_summary(local_opt):
    """Simplifiy option dict creation"""
    d = OrderedDict()
    d["Package name"] = pkg_name
    d["Package version"] = pkg_version
    d["Timestamp"] = str(datetime.datetime.now())
    for i, j in local_opt.items():
        d[i] = j
    return d

def assert_column(df, column_name, df_name=''):
        assert column_name in df.columns, f"'{column_name}' column is missing in {df_name}"
        assert not df[column_name].isnull().all(), f"'{column_name}' column is empty in {df_name}"


def call_bash_script(bash_script_path, *args):
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the path to the bash script
    bash_script_path = os.path.join(script_dir, bash_script_path)

    # Change the permissions of the bash script to make it executable
    os.chmod(bash_script_path, stat.S_IRWXU)

    # Create the command by joining the script path and arguments
    command = [bash_script_path] + list(args)

    # Call the command using subprocess.Popen
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Get the output and error messages
    stdout, stderr = process.communicate()

    # Decode the output and error messages
    stdout = stdout.decode('utf-8')
    stderr = stderr.decode('utf-8')

    # Return the output and error messages
    return stdout, stderr


def str_join(l, sep="\t", line_end=""):
    """Join a list of mixed types into a single str"""
    s = sep.join(map(str, l)) + line_end
    return s


def list_to_str(l):
    """Generate a string from any list"""
    return str(json.dumps(l)).replace(" ", "")


def str_to_list(s, parse_int=None, parse_float=None):
    """Generate a list from a string"""
    return json.loads(s, parse_int=parse_int, parse_float=parse_float)


def all_in(l1, l2):
    """Check if all element in l1 are in l2"""
    s1 = set(l1)
    s2 = set(l2)
    if s1.difference(s2):
        return False
    else:
        return True


def file_readable(fn, **kwargs):
    """Check if the file is readable"""
    return os.path.isfile(fn) and os.access(fn, os.R_OK)


def dir_writable(fn, **kwargs):
    """Check if the file is readable"""
    if not os.path.isdir(fn):
        fn = os.path.dirname(fn)
    return os.path.dirname(fn) and os.access(fn, os.W_OK)


def mkdir(fn, exist_ok=False):
    """ Create directory recursivelly. Raise IO error if path exist or if error at creation """
    try:
        os.makedirs(fn, exist_ok=exist_ok)
    except:
        raise ironsightError("Error creating output folder `{}`".format(fn))


def mkbasedir(fn, exist_ok=False):
    """ Create directory for a given file recursivelly. Raise IO error if path exist or if error at creation """
    dir_fn = os.path.dirname(fn)
    if dir_fn:
        mkdir(dir_fn, exist_ok=True)


def dict_to_str(d, sep="\t", nsep=0, exclude_list=[]):
    """ Transform a multilevel dict to a tabulated str """
    m = ""
    
    if isinstance(d, Counter):
        for i, j in d.most_common():
            if not i in exclude_list:
                m += "{}{}: {:,}\n".format(sep * nsep, i, j)
    
    else:
        for i, j in d.items():
            if not i in exclude_list:
                if isinstance(j, dict):
                    j = dict_to_str(j, sep=sep, nsep=nsep + 1)
                    m += "{}{}\n{}".format(sep * nsep, i, j)
                else:
                    m += "{}{}: {}\n".format(sep * nsep, i, j)
    if not m:
        return ""
    else:
        return m[:-1]


def iter_idx_tuples(df):
    for idx, line in zip(df.index, df.itertuples(index=False, name="line")):
        yield (idx, line)


def doc_func(func):
    """Parse the function description string"""
    
    if inspect.isclass(func):
        func = func.__init__
    
    docstr_list = []
    for l in inspect.getdoc(func).split("\n"):
        l = l.strip()
        if l:
            if l.startswith("*"):
                break
            else:
                docstr_list.append(l)
    
    return " ".join(docstr_list)


def make_arg_dict(func):
    """Parse the arguments default value, type and doc"""
    
    # Init method for classes
    if inspect.isclass(func):
        func = func.__init__
    
    if inspect.isfunction(func) or inspect.ismethod(func):
        # Parse arguments default values and annotations
        d = OrderedDict()
        for name, p in inspect.signature(func).parameters.items():
            if not p.name in ["self", "cls"]:  # Object stuff. Does not make sense to include in doc
                d[name] = OrderedDict()
                if not name in ["kwargs", "args"]:  # Include but skip default required and type
                    # Get Annotation
                    if p.annotation != inspect._empty:
                        d[name]["type"] = p.annotation
                    # Get default value if available
                    if p.default == inspect._empty:
                        d[name]["required"] = True
                    else:
                        d[name]["default"] = p.default
        
        # Parse the docstring in a dict
        docstr_dict = OrderedDict()
        lab = None
        for l in inspect.getdoc(func).split("\n"):
            l = l.strip()
            if l:
                if l.startswith("*"):
                    lab = l[1:].strip()
                    docstr_dict[lab] = []
                elif lab:
                    docstr_dict[lab].append(l)
        
        # Concatenate and copy doc in main dict
        for name in d.keys():
            if name in docstr_dict:
                d[name]["help"] = " ".join(docstr_dict[name])
        return d


def arg_from_docstr(parser, func, arg_name, short_name=None):
    """Get options corresponding to argument name from docstring and deal with special cases"""
    
    if short_name:
        arg_names = ["-{}".format(short_name), "--{}".format(arg_name)]
    else:
        arg_names = ["--{}".format(arg_name)]
    
    arg_dict = make_arg_dict(func)[arg_name]
    if "help" in arg_dict:
        if "default" in arg_dict:
            if arg_dict["default"] == "" or arg_dict["default"] == []:
                arg_dict["help"] += " (default: None)"
            else:
                arg_dict["help"] += " (default: %(default)s)"
        else:
            arg_dict["help"] += " (required)"
        
        if "type" in arg_dict:
            arg_dict["help"] += " [%(type)s]"
    
    # Special case for boolean args
    if arg_dict["type"] == bool:
        if arg_dict["default"] == False:
            arg_dict["action"] = "store_true"
            del arg_dict["type"]
        elif arg_dict["default"] == True:
            arg_dict["action"] = "store_false"
            del arg_dict["type"]
    
    # Special case for lists args
    elif isinstance(arg_dict["type"], list):
        arg_dict["nargs"] = "*"
        arg_dict["type"] = arg_dict["type"][0]
    
    parser.add_argument(*arg_names, **arg_dict)


def head(fp, n=10, sep="\t", max_char_col=50, comment=None):
    """
    Emulate linux head cmd. Handle gziped files and bam files
    * fp
        Path to the file to be parse.
    * n
        Number of lines to print starting from the begining of the file (Default 10)
    """
    line_list = []
    
    # Get lines
    try:
        open_fun, open_mode = (gzip.open, "rt") if fp.endswith(".gz") else (open, "r")
        with open_fun(fp, open_mode) as fh:
            line_num = 0
            while line_num < n:
                l = next(fh).strip()
                if comment and l.startswith(comment):
                    continue
                if sep:
                    line_list.append(l.split(sep))
                else:
                    line_list.append(l)
                line_num += 1
    
    except StopIteration:
        pass
    
    # Add padding if sep given
    if sep:
        try:
            # Find longest elem per col
            col_len_list = [0 for _ in range(len(line_list[0]))]
            for ls in line_list:
                for i in range(len(ls)):
                    len_col = len(ls[i])
                    if len_col > max_char_col:
                        col_len_list[i] = max_char_col
                    elif len_col > col_len_list[i]:
                        col_len_list[i] = len_col
            
            # Add padding
            line_list_tab = []
            for ls in line_list:
                s = ""
                for i in range(len(ls)):
                    len_col = col_len_list[i]
                    len_cur_col = len(ls[i])
                    if len_cur_col <= len_col:
                        s += ls[i] + " " * (len_col - len_cur_col) + " "
                    else:
                        s += ls[i][0 : len_col - 3] + "..."
                line_list_tab.append(s)
            line_list = line_list_tab
        
        # Fall back to non tabulated display
        except IndexError:
            return head(fp=fp, n=n, sep=None)
    
    for l in line_list:
        print(l)
    print()


def stdout_print(*args):
    """
    Emulate print but uses sys stdout instead. It could sometimes be useful in specific situations where print
    is in is not behaving optimaly (like with tqdm for example)
    """
    s = " ".join([str(i) for i in args])
    sys.stdout.write(s)
    sys.stdout.flush()


def get_logger(name=None, verbose=False, quiet=False):
    """Multilevel colored log using colorlog"""
    
    # Define conditional color formatter
    formatter = colorlog.LevelFormatter(
        fmt={
            "DEBUG": "%(log_color)s\t[DEBUG]: %(msg)s",
            "INFO": "%(log_color)s\t%(msg)s",
            "WARNING": "%(log_color)s## %(msg)s ##",
            "ERROR": "%(log_color)sERROR: %(msg)s",
            "CRITICAL": "%(log_color)sCRITICAL: %(msg)s",
        },
        log_colors={
            "DEBUG": "white",
            "INFO": "green",
            "WARNING": "bold_blue",
            "ERROR": "bold_red",
            "CRITICAL": "bold_purple",
        },
        reset=True,
    )
    
    # Define logger with custom formatter
    logging.basicConfig(format="%(message)s")
    logging.getLogger().handlers[0].setFormatter(formatter)
    log = logging.getLogger(name)
    
    # Define logging level depending on verbosity
    if verbose:
        log.setLevel(logging.DEBUG)
    elif quiet:
        log.setLevel(logging.WARNING)
    else:
        log.setLevel(logging.INFO)
    
    return log


def log_dict(d, logger, header="", indent="\t", level=1):
    """ log a multilevel dict """
    if header:
        logger(header)
    if isinstance(d, Counter):
        for i, j in d.most_common():
            logger("{}{}: {:,}".format(indent * level, i, j))
    else:
        for i, j in d.items():
            if isinstance(j, dict):
                logger("{}{}".format(indent * level, i, j))
                log_dict(j, logger, level=level + 1)
            else:
                logger("{}{}: {}".format(indent * level, i, j))


def log_list(l, logger, header="", indent="\t"):
    """ log a list """
    if header:
        logger(header)
    for i in l:
        logger("{}*{}".format(indent, i))


# ~~~~~~~~~~~~~~CUSTOM EXCEPTION AND WARN CLASSES~~~~~~~~~~~~~~#
class ironsightError(Exception):
    """ Basic exception class for ironsight package """
    
    pass

