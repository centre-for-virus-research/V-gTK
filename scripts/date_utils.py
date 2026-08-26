from datetime import datetime

from dateutil import parser
import re

# ISO 8601, whole or partial: 2020-03-15, 2020-03, and the sloppy 2020-3-5.
# Handled explicitly and BEFORE the dateutil path below, because the heuristics
# there cannot see a numeric month.
_ISO_RE = re.compile(r'^(\d{4})-(\d{1,2})(?:-(\d{1,2}))?$')


def split_date_components(date_str):
    """Split a GenBank collection_date into day / month / year components.

    Missing components come back as ``''`` rather than a guess, so that a
    year-only date never becomes 1 January.

    ISO 8601 is handled explicitly. The generic path below infers which
    components are present with three regexes, and its month test only matches
    an alphabetic month name (Jan, Feb, ...). That is fine for GenBank's legacy
    ``15-Mar-2020`` but silently dropped the day AND month of every ISO date:
    ``2020-03-15`` came back as year-only because '03' is not 'Mar', and because
    the day test inspects the text before the first '-', which for an ISO date is
    the year. INSDC is standardising on ISO, so the most correct format was the
    one losing the most precision - 5,598 rows in a single HCV database had a
    full date in the source and a blank day and month in the built database.
    """
    try:
        date_str_clean = date_str.strip()
        if not date_str_clean:
            return {'day': '', 'month': '', 'year': ''}

        iso = _ISO_RE.match(date_str_clean)
        if iso:
            year, month, day = iso.group(1), iso.group(2), iso.group(3)
            try:
                # datetime validates the calendar for us: 2020-02-30 and
                # 2020-13-01 raise rather than silently rolling over.
                datetime(int(year), int(month), int(day) if day else 1)
            except ValueError:
                # The year is still known and still usable; only the parts that
                # failed to validate are dropped.
                return {'day': '', 'month': '', 'year': int(year)}
            return {
                'day': int(day) if day else '',
                'month': int(month),
                'year': int(year),
            }

        mmddyy_pattern = r'^\d{1,2}/\d{1,2}/\d{2}$'
        if re.match(mmddyy_pattern, date_str_clean):
            dt = parser.parse(date_str_clean, dayfirst=False, yearfirst=False)
            return {'day': dt.day, 'month': dt.month, 'year': dt.year}

        dt = parser.parse(date_str_clean, default=parser.parse("01-Jan-1900"))

        # Infer presence based on patterns
        has_day = bool(re.search(r'\b\d{1,2}\b', date_str_clean.split('-')[0]))
        has_month = bool(re.search(r'\b(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\b', date_str_clean, re.IGNORECASE))
        has_year = bool(re.search(r'\b\d{4}\b', date_str_clean))

        # A day without a month, or a day/month without a year, is not a partial
        # date - it is a misparse. '8/25/2022' used to yield day=25 with no month
        # (the month test only knows names), and '15-Mar-20' yielded day and
        # month with no year at all. Neither is safe to store: downstream code
        # cannot tell a real partial date from a mangled one.
        if not has_year:
            return {'day': '', 'month': '', 'year': ''}
        if has_day and not has_month:
            return {'day': '', 'month': '', 'year': dt.year}

        return {
            'day': dt.day if has_day else '',
            'month': dt.month if has_month else '',
            'year': dt.year if has_year else ''
        }

    except Exception:
        return {'day': '', 'month': '', 'year': ''}
