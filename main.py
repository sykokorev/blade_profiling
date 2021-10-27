# -*- coding: utf-8 -*-
import subprocess
import sys
import os


if __name__ == '__main__':
    root_dir = os.getcwd()
    autobalde = os.path.join(root_dir, 'main_autoblade.py')
    p = subprocess.Popen([sys.executable, autobalde])
    p.wait()
    nx_script = os.path.join(root_dir, 'main_nx.py')
    nx = r'C:\Program Files\Siemens\NX 12.0\NXBIN\run_journal.exe '
    nx += nx_script
    p = subprocess.call(nx)
