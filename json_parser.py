'''Parse data in JSON format.'''

import sys
import json

for line in open(sys.argv[1]):
    jsn = json.loads(line.strip())
    print '%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s' % \
                                        (jsn['votes']['cool'],
                                        jsn['votes']['funny'],
                                        jsn['votes']['useful'],
                                        jsn['user_id'],
                                        jsn['name'],
                                        jsn['url'],
                                        jsn['average_stars'],
                                        jsn['review_count'],
                                    )
