'''
Simple queue listing tool.

Prints such that the bottom task is that most expected to finish
(highest priority in queue, or earliest submitted).

For each task, the
ID, number of cores, status, and date started/submitted
are shown, followed by the task name.

Use --unqueued alongside a list of task names to maintain your queue:

qstat --unqueued --list path/to/names.txt | awk '{print "qsub "$0".sh"}'
'''

from argparse import ArgumentParser
import os
from sys import argv
from time import sleep, time
import xml.etree.ElementTree as ET

from colorama import Fore, Style


class Queue:
    '''
    State-machine for qstat to recall tasks not queued or running
    (and tasks that are complete, if ran on a loop)
    '''

    def __init__(self, namepath):
        self.lastseen = dict()
        self.finished = []

        self.names = []
        if namepath is not None and os.path.isfile(namepath):
            with open(namepath) as f:
                for l in f.readlines():
                    l = l.strip()
                    if l and not l.startswith('#'):
                        self.names.append(l)

    def jobs(self):
        data = os.popen('qstat -xml').read()
        tree = ET.XML(data)

        queued = list(tree.find('job_info'))
        queued.sort(key=lambda j: float(j.find('JAT_prio').text))

        running = list(tree.find('queue_info'))
        running.sort(key=lambda j: (
            j.find('JAT_start_time').text,
            j.find('JB_job_number').text
        ), reverse=True)

        return queued + running

    def run(self):
        jobs = self.jobs()

        # to maintain self.finished
        unseen = dict(self.lastseen)
        nowseen = dict()
        # to find tasks not in self.names
        unseen_names = set(self.names or set())

        for job in jobs:
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

            # handle self.finished and self.names
            if state != 'dr':
                # if marked for deletion, don't put into "seen tasks"
                nowseen[ID] = (name, text)
            if ID in unseen:
                unseen.pop(ID)
            if name in unseen_names:
                unseen_names.remove(name)

        self.lastseen = nowseen

        for ID, (name, text) in unseen.items():
            if name in unseen_names:
                unseen_names.remove(name)
            self.finished.insert(0, text)

        for job in self.finished:
            print(Fore.GREEN + job.replace(' r  ', ' d  '))
        for name in unseen_names:
            print(
                Fore.LIGHTMAGENTA_EX,
                ' ' * 18,
                'not in queue or running:',
                Style.BRIGHT + name + Style.RESET_ALL
            )

    def unseen_names(self):
        desired = set(self.names)
        seen = set(job.find('JB_name').text for job in self.jobs())
        return desired - seen


parser = ArgumentParser(description=__doc__)
parser.add_argument('--loop', action='store_true',
                    help='loop every 10 seconds')
parser.add_argument('--list', metavar='path', type=str, default=None,
                    help='''
path to a list of task names that should be in queue (are also printed)
''')
parser.add_argument('--unqueued', action='store_true',
                    help='if --list set, print only names of tasks NOT in queue')

parser.add_argument('--ids', action='store_true',
                    help='list IDs separated by spaces')

if __name__ == '__main__':
    args = parser.parse_args()
    queue = Queue(args.list)
    if args.ids:
        print(' '.join([
            job.find('JB_job_number').text
            for job in queue.jobs()
        ]))
    elif args.unqueued:
        for name in queue.unseen_names():
            print(name)
    elif args.loop:
        while True:
            try:
                os.system('cls' if os.name == 'nt' else 'clear')
                queue.run()
                sleep(10)
            except KeyboardInterrupt:
                break
    else:
        queue.run()
