#!/usr/bin/python
#
import sys
import time
from marzhauser import API



def main():
    '''
    load the API and run a couple of functions
    '''
    mh_api = API()
    MX = float(sys.argv[1])

    if mh_api.getStatus():
        print('The stage is functioning well.')
        current_pos = mh_api.position()
        print(current_pos)
        mh_api.setVelocity('x', 0.5 )
        mh_api.goRelative(MX, 0.0)
        mh_api.setVelocity('x', 5.0)

        mh_api.shutDown()

    else:
        print("Failed to initialize the stage.")


if __name__ == '__main__':
    main()
