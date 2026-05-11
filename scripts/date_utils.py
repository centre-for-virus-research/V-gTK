from dateutil import parser
import re

def split_date_components(date_str):
    try:
        date_str_clean = date_str.strip()
        if not date_str_clean:
            return {'day': '', 'month': '', 'year': ''}


        mmddyy_pattern = r'^\d{1,2}/\d{1,2}/\d{2}$'
        if re.match(mmddyy_pattern, date_str_clean):
            dt = parser.parse(date_str_clean, dayfirst=False, yearfirst=False)
            return {'day': dt.day,'month': dt.month,'year': dt.year}

        dt = parser.parse(date_str_clean, default=parser.parse("01-Jan-1900"))

        # Infer presence based on patterns
        has_day = bool(re.search(r'\b\d{1,2}\b', date_str_clean.split('-')[0]))
        has_month = bool(re.search(r'\b(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\b', date_str_clean, re.IGNORECASE))
        has_year = bool(re.search(r'\b\d{4}\b', date_str_clean))

        return {
            'day': dt.day if has_day else '',
            'month': dt.month if has_month else '',
            'year': dt.year if has_year else ''
        }

    except Exception:
        return {'day': '', 'month': '', 'year': ''}

#print(split_date_components("14-Sep-2010"))  # {'day': 14, 'month': 9, 'year': 2010}
#print(split_date_components("2013"))         # {'day': '', 'month': '', 'year': 2013}
#print(split_date_components("Aug-2013"))     # {'day': '', 'month': 8, 'year': 2013}
#print(split_date_components("8/25/2022"))     # {'day': '', 'month': '', 'year': ''}
