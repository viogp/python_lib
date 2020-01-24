import calendar,datetime

def year_days(year):
    d1 = datetime.date(year, 1, 1) 
    d2 = datetime.date(year + 1, 1, 1) 
    return (d2 - d1).days 

def find_month_day(year,yday):
    this_year = year_days(year)
    if ((yday<1) or (yday>this_year)):
        print('STOP: Day {} is outside year {}'.format(yday,year))
        sys.exit()

    last_month = 0 ; adday = 0 ; i = 1
    month_found = False
    while (not month_found):
        # Returns weekday of first day of the month and 
        # number of days in month
        adday += calendar.monthrange(year, i)[1]
        if (yday<=adday):
            month = i            
            day = yday - last_month
            month_found = True
        else:
            last_month = adday 
            i += 1 

    return month, day

def find_week(year,yday):
    this_year = year_days(year)
    if ((yday<1) or (yday>this_year)):
        print('STOP: Day {} is outside year {}'.format(yday,year))
        sys.exit()

    month, day = find_month_day(year,yday)
    # isocalendar: year, week number, day number
    week = datetime.date(year, month, day).isocalendar()[1]
    yweek = datetime.date(year, month, day).isocalendar()[0]
    #if (yweek < year):
    #    print('WARNING: week {} corresponds to the previous year {}'.format(week,yweek))

    return week, yweek

if __name__ == '__main__':
    year = int(input('Enter a year: '))
    days = year_days(year)
    print('Year {} had {} days'.format(year,days))

    inday = int(input('Enter a day (1-365): '))
    month, day = find_month_day(year,inday)
    print('Day {} = {}/{}/{}'.format(inday,day,month,year))    
