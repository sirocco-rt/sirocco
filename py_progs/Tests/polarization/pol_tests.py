#!/usr/bin/env python
'''
Run the polarization tests in py_progs/Tests/polarization and summarize whether
they pass.


Command line usage (if any):

    usage::

        pol_tests.py [-h] [-full] [-np N] [-exe SIROCCO] [-w2t WINDSAVE2TABLE]
                     [-tests TEST [TEST ...]] [-check_only] workdir

    workdir     the directory in which to run the tests; it is created if
                necessary, and must not be inside $SIROCCO
    -full       run the full tests (several hours on 8 processes) rather
                than the quick ones (about 10 minutes)
    -np         number of MPI processes (default 8; 1 runs without mpirun)
    -exe        the Sirocco executable to test (default sirocco)
    -w2t        a windsave2table built from the same code as -exe (default
                windsave2table); needed by the cvcone test
    -tests      run only these tests (default: brown_mclean hillier
                sr_shell cvcone)
    -check_only do not run Sirocco; repeat the analysis and the checks on
                the results already in workdir

Description:

    Each test is a subdirectory of py_progs/Tests/polarization.  It holds scripts
    that write the models, an independent prediction, and testdef.py, which
    says which models to run in the quick and the full tiers, which analysis
    scripts to run afterwards, and what to check.

    For each test, the files are copied to workdir/test, Setup_Sirocco_Dir
    is run there, the models are made and run with -d, and then the checks
    are made.  Nothing is written in the repository.  The results go to
    workdir/Summary.txt (and the screen): for each check, the value, the
    range it must lie in, and PASS or FAIL, plus the executable and its
    version, so that runs of different versions can be compared.  The
    model outputs and the comparison plots stay in workdir/test.

Primary routines:

    doit    runs the tests; pol_tests.doit(workdir, ...) can be called
            from a script

Notes:

    $SIROCCO must be set.  Requires numpy and matplotlib.

History:

    261018 ksl Coding begun, with Claude Code

'''
import datetime
import importlib.util
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = ['brown_mclean', 'hillier', 'sr_shell', 'cvcone']
SHARED = ['h_only_master.dat', 'h_elem_ions.dat']


def version(exe, cwd='.'):
    '''The version, git hash and any uncommitted changes from exe --version, to record what was tested'''
    try:
        proc = subprocess.run([exe, '--version'], capture_output=True, text=True, timeout=60, cwd=cwd)
        lines = [one.strip() for one in proc.stdout.splitlines()]
        keep = [one for one in lines if 'Version' in one or 'hash' in one or 'uncommitted' in one]
        return keep if keep else ['no version information from %s --version' % exe]
    except (OSError, subprocess.TimeoutExpired) as err:
        return ['could not run %s --version: %s' % (exe, err)]


def setup(test, testdir):
    '''Copy a test, and the shared atomic data, to testdir and link the Sirocco data there'''
    os.makedirs(testdir, exist_ok=True)
    src = os.path.join(HERE, test)
    for one in os.listdir(src):
        if one != '__pycache__' and os.path.isfile(os.path.join(src, one)):
            shutil.copy(os.path.join(src, one), testdir)
    for one in SHARED:
        shutil.copy(os.path.join(HERE, one), testdir)
    if not os.path.exists(os.path.join(testdir, 'data')):
        subprocess.run([os.path.join(os.environ['SIROCCO'], 'bin', 'Setup_Sirocco_Dir')], cwd=testdir,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def load(test, testdir):
    '''Import the testdef.py of a test from its working directory'''
    spec = importlib.util.spec_from_file_location('testdef_' + test, os.path.join(testdir, 'testdef.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def completed(root, testdir):
    '''True if the .sig file says that the run of root completed'''
    try:
        with open(os.path.join(testdir, root + '.sig')) as f:
            lines = f.read().splitlines()
        return len(lines) > 0 and 'COMPLETE' in lines[-1]
    except IOError:
        return False


def run_models(tdef, tier, testdir, exe, np, w2t, log):
    '''Make and run the models of one test; returns a list of problems'''
    problems = []
    for model in tdef.models(tier):
        root = model['root']
        if model.get('make'):
            proc = subprocess.run([sys.executable] + model['make'], cwd=testdir, capture_output=True, text=True)
            if proc.returncode:
                problems.append('%s: %s failed: %s' % (root, model['make'][0], proc.stderr.strip()[-300:]))
                continue
        command = [exe, '-d'] + model.get('switches', []) + [root]
        if np > 1:
            command = ['mpirun', '-np', str(np)] + command
        log('    running %s' % ' '.join(command))
        with open(os.path.join(testdir, root + '.stdout'), 'w') as out:
            subprocess.run(command, cwd=testdir, stdout=out, stderr=subprocess.STDOUT)
        if not completed(root, testdir):
            problems.append('%s: the run did not complete (see %s.stdout and %s.sig)' % (root, root, root))
            continue
        if model.get('w2t'):
            proc = subprocess.run([w2t, root], cwd=testdir, capture_output=True, text=True)
            if proc.returncode or not os.path.exists(os.path.join(testdir, root + '.master.txt')):
                problems.append('%s: %s could not read the wind_save file; it must be built from the same code '
                                'as the Sirocco being tested (use -w2t)' % (root, w2t))
    return problems


def run_post(tdef, tier, testdir, log):
    '''Run the analysis scripts of one test; returns a list of problems'''
    problems = []
    for command in tdef.post(tier):
        log('    %s' % ' '.join([os.path.basename(command[0])] + command[1:]))
        logname = os.path.splitext(os.path.basename(command[0]))[0] + '.log'
        with open(os.path.join(testdir, logname), 'a') as out:
            proc = subprocess.run([sys.executable] + command, cwd=testdir, stdout=out, stderr=subprocess.STDOUT)
        if proc.returncode:
            problems.append('%s failed (see %s)' % (' '.join(command), logname))
    return problems


def run_checks(tdef, tier, testdir):
    '''Make the checks of one test in its directory; returns (results, problems)'''
    cwd = os.getcwd()
    sys.path.insert(0, testdir)
    try:
        os.chdir(testdir)
        return tdef.checks(tier), []
    except Exception as err:
        return [], ['the checks could not be made: %s: %s' % (type(err).__name__, err)]
    finally:
        os.chdir(cwd)
        sys.path.remove(testdir)


def doit(workdir, tier='quick', np=8, exe='sirocco', w2t='windsave2table', tests=TESTS, check_only=False):
    '''
    Run the tests in workdir and write workdir/Summary.txt.  Returns True if
    every check passed.
    '''
    workdir = os.path.abspath(workdir)
    os.makedirs(workdir, exist_ok=True)
    lines = []

    def log(line):
        print(line)
        lines.append(line)

    log('Polarization tests (%s), %s' % (tier, datetime.datetime.now().strftime('%Y-%m-%d %H:%M')))
    log('Executable: %s, %d process%s; windsave2table: %s' % (exe, np, 'es' if np > 1 else '', w2t))
    for one in version(exe, workdir):
        log('    %s' % one)
    if check_only:
        log('(-check_only: Sirocco was not run; checks of the results already in %s)' % workdir)
    log('')

    npass = nfail = 0
    for test in tests:
        testdir = os.path.join(workdir, test)
        t0 = time.time()
        if not check_only:
            setup(test, testdir)
        tdef = load(test, testdir)
        log('%s: %s' % (test, tdef.DESCRIPTION))
        problems = []
        if not check_only:
            problems += run_models(tdef, tier, testdir, exe, np, w2t, log)
        if not problems:
            problems += run_post(tdef, tier, testdir, log)
        results, more = run_checks(tdef, tier, testdir) if not problems else ([], [])
        problems += more
        for label, value, lo, hi in results:
            if lo is None:
                log('    INFO  %-78s %+9.4f' % (label, value))
            else:
                ok = lo <= value <= hi
                npass += ok
                nfail += not ok
                log('    %s  %-78s %+9.4f   [%g, %g]' % ('PASS' if ok else 'FAIL', label, value, lo, hi))
        for one in problems:
            nfail += 1
            log('    FAIL  %s' % one)
        log('    (%.1f min; results in %s)' % ((time.time() - t0) / 60, testdir))
        log('')

    log('%d checks passed, %d failed: %s' % (npass, nfail, 'ALL PASSED' if nfail == 0 else 'SOME FAILED'))
    with open(os.path.join(workdir, 'Summary.txt'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print('Summary written to %s' % os.path.join(workdir, 'Summary.txt'))
    return nfail == 0


def steer(argv):
    '''
    Process the command line; returns the exit status (0 if every check passed)
    '''

    workdir = ''
    tier = 'quick'
    nproc = 8
    exe = 'sirocco'
    w2t = 'windsave2table'
    tests = TESTS
    check_only = False

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return 0
        elif argv[i] == '-full':
            tier = 'full'
        elif argv[i] == '-np':
            i += 1
            nproc = int(argv[i])
        elif argv[i] == '-exe':
            i += 1
            exe = argv[i]
        elif argv[i] == '-w2t':
            i += 1
            w2t = argv[i]
        elif argv[i] == '-tests':
            tests = []
            while i + 1 < len(argv) and argv[i + 1] in TESTS:
                i += 1
                tests.append(argv[i])
            if len(tests) == 0:
                print('Error: -tests needs one or more of %s' % TESTS)
                return 1
        elif argv[i] == '-check_only':
            check_only = True
        elif argv[i][0] == '-':
            print('Error: Unknown switch ---  %s' % argv[i])
            return 1
        else:
            workdir = argv[i]
        i += 1

    if workdir == '':
        print(__doc__)
        return 1
    if 'SIROCCO' not in os.environ:
        print('Error: $SIROCCO must be set')
        return 1
    repo = os.path.realpath(os.environ['SIROCCO'])
    if os.path.realpath(os.path.abspath(workdir)).startswith(repo + os.sep):
        print('Error: run the tests outside the repository, not in %s' % repo)
        return 1
    ok = doit(workdir, tier, nproc, exe, w2t, tests, check_only)
    return 0 if ok else 1


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        sys.exit(steer(sys.argv))
    else:
        print(__doc__)
