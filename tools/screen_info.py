
import sys
import time
from time import sleep


# prints a phrase which get continuously
# updated on the same line
def update_print(phrase, appendix=' ', percent=False):
    sys.stdout.write('\r')
    if percent==True:
        sys.stdout.write(phrase + '%d%%' % appendix )
        if appendix == 100:
            print '\n'
    if percent==False:
        sys.stdout.write(phrase + appendix )
        if appendix == 'done':
            print '\n'
    sys.stdout.flush()
    sleep(1.e-6)
