'''
Simple queue listing tool.

Prints such that the bottom task is that most expected to finish
(highest priority in queue, or earliest submitted).

For each task, the
ID, number of cores, status, and date started/submitted
are shown, followed by the task name.
If --names or --ids is provided, only those are given.

Use --unqueued, --queued, --running (or some of -uqr) to list only certain jobs.
By default, all are listed.

Examples of use:

qstat -un --list path/to/names.txt | awk '{print "qsub "$0".sh"}'
qdel $(qstat -qi)
'''

from argparse import ArgumentParser
import os
from sys import argv
from time import sleep, time
import xml.etree.ElementTree as ET

from colorama import Fore, Style

def run(args):
    names = []
    if args.list is not None and os.path.isfile(args.list):
        with open(args.list) as f:
            for l in f.readlines():
                l = l.strip()
                if l and not l.startswith('#'):
                    names.append(l)
    #
    # Parse...
    #
    tree = ET.XML(os.popen('qstat -xml').read())

    # Get list of jobs
    jobs_show = []
    jobs_seen = []

    use_default = not (args.unqueued or args.queued or args.running)
    
    queued = list(tree.find('job_info'))
    queued.sort(key=lambda j: float(j.find('JAT_prio').text))
    jobs_seen += queued
    if use_default or args.queued:
        jobs_show += queued

    running = list(tree.find('queue_info'))
    running.sort(key=lambda j: (
        j.find('JAT_start_time').text,
        j.find('JB_job_number').text
    ), reverse=True)
    jobs_seen += running
    if use_default or args.running:
        jobs_show += running
    
    # find tasks not in --list
    unseen_names = set()
    if use_default or args.unqueued:
        unseen_names = set(names or set())
        for job in jobs_seen:
            name = job.find('JB_name').text
            if name in unseen_names:
                unseen_names.remove(name)

    #
    # ...and display
    #
    
    if args.ids:
        print(' '.join([
            job.find('JB_job_number').text
            for job in jobs_show
        ]))
    elif args.names:
        for job in jobs_show:
            print(job.find('JB_name').text)
        for name in unseen_names:
            print(name)
    else:
        for job in jobs_show:
            def F(length, *names):
                for x in names:
                    n = job.find(x)
                    if n is not None:
                        break
                else:
                    return ' ' * length
                return n.text.ljust(length)
            state = F(3, 'state')
            color = Fore.YELLOW
            verb = 'submitted '
            if state.strip() == 'dr':
                verb = '  started '
                color = Fore.RED
            elif 'q' in state:
                color = Fore.RED
            elif 'r' in state:
                color = Fore.CYAN
                verb = '  started '

            ID = F(7, 'JB_job_number')
            slots = F(2, 'slots')
            stateS = Style.BRIGHT + state + Style.RESET_ALL + color
            date = F(19, 'JAT_start_time', 'JB_submission_time')
            date = Style.DIM + verb + date.replace('T', ' ')

            name = F(0, 'JB_name')
            nameS = Style.RESET_ALL + color + Style.BRIGHT
            nameS += name + Style.RESET_ALL

            text = f"{ID} {slots} {stateS} {date} {nameS}"
            print(color + text)
        
        for name in unseen_names:
            print(
                Fore.LIGHTMAGENTA_EX,
                ' ' * 19,
                'not in queue or running',
                Style.BRIGHT + name + Style.RESET_ALL
            )


parser = ArgumentParser(description=__doc__)
parser.add_argument('--list', metavar='path', type=str, default=None,
                    help='''
path to a list of task names that should be in queue (are also printed)
''')

# todo: more dynamic args (eg --unqueued becomes -un aka --unqueued --names)

parser.add_argument('--unqueued', '-u', action='store_true',
                    help='only show names in --list but not queued or running')
parser.add_argument('--queued', '-q', action='store_true',
                    help='only show tasks queued to run')
parser.add_argument('--running', '-r', action='store_true',
                    help='only show tasks running')

parser.add_argument('--ids', '-i', action='store_true',
                    help='list job IDs separated by spaces')
parser.add_argument('--names', '-n', action='store_true',
                    help='list names of jobs or potential tasks on each line')

if __name__ == '__main__':
    args = parser.parse_args()
    run(args)
