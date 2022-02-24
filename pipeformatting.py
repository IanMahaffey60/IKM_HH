import csv

with open(r"L:\LAProj\10067.006 - Building D-2 - Victory Logistics\Civil\Hydrology\calculations\hydro\results\pipe_report.csv", 'r') as f:
    filer = csv.reader(f, delimiter=',')
    with open(r"L:\LAProj\10067.006 - Building D-2 - Victory Logistics\Civil\Hydrology\calculations\hydro\results\report_formatting1.txt", 'w') as f:
        for num, row in enumerate(filer):
            if num == 0:
                pass
            else:
                print (f'Name = {row[0]}', file=f)
                print (f'Q(5)={round(float(row[4]),1)} cfs ~ Q(25)={round(float(row[5]),1)} cfs', file=f)
                print (f'V(5)={round(float(row[6]),1)} fps ~ V(25)={round(float(row[7]),1)} fps', file=f)
                print (f'Q(cap)={round(float(row[3]),1)} cfs\n', file=f)
