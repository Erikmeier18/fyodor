import astropy.units as u
from astropy.coordinates import EarthLocation

__all__ = ['Location', 'cerro_paranal', 'san_pedro_martir']

cp = dict(
    latdeg = -24.6230,
    londeg = -70.4025,
    site = 'Cerro Paranal'
)

spm = dict(
    latdeg = 30.9058267,
    londeg = -115.4254954,
    site = 'San Pedro Mártir'
)


class Location(object):
    def __init__(self, latdeg, londeg, site):
        self.latdeg = latdeg
        self.londeg = londeg
        self.site = site

    def to_EarthLocation(self):
        return EarthLocation.from_geodetic(lon=self.londeg * u.deg,
                                           lat=self.latdeg * u.deg)


cerro_paranal = Location(**cp)
san_pedro_martir = Location(**spm)
