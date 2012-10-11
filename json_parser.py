'''Parse data in JSON format.'''

import sys
import csv
import json

with open('output.csv', 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    for line in open(sys.argv[1]):
        jsn = json.loads(line.strip())
        writer.writerow([jsn['votes']['cool'],
                            jsn['votes']['funny'],
                            jsn['votes']['useful'],
                            jsn['user_id'],
                            jsn['name'],
                            jsn['url'],
                            jsn['average_stars'],
                            jsn['review_count'],
                        ])
