from sunpy.net import Fido, attrs as a
import astropy.units as u
import os
import datetime

# Select the year and by changing the year from 2011 to 2023, 
#then it is able to download the solar image for that time period.

start_year = 2012
end_year = 2012


for year in range(start_year, end_year + 1):
    for month in range(1, 13): 
        # 
        if month == 12:
            start_date = datetime.datetime(year, month, 1)
            end_date = datetime.datetime(year, month, 31)
        else:
            
            import calendar
            _, last_day = calendar.monthrange(year, month)
            start_date = datetime.datetime(year, month, 1)
            end_date = datetime.datetime(year, month, last_day)

        # Setting up the conditions, selecting the intensity image of BBSO
        timerange = a.Time(start_date, end_date)
        results = Fido.search(timerange,
                              a.Source.gong,
                              a.Provider.nso,
                              a.Physobs.intensity,
                              a.Instrument.bigbear,
                              a.Sample(1 * u.day))
                              
        ######setting the dictionary
        directory = "/file/HAtest2/{year}/{month:02d}".format(year=year, month=month)  
        os.makedirs(directory, exist_ok=True)

        #acquire the data
        files = Fido.fetch(results, path=os.path.join(directory, "{instrument}_{file}"))
        
        print(f'Downloaded {len(files)} files for {year}-{month:02d}')
