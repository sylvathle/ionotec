import scipy.constants as csts

# Channels of the glonass constellation per satellite ("k" values)
 #https://glonass-iac.ru/en/cus/

gps_f1, gps_f2, gps_f5 = 1575.42 * 1e6, 1227.60 * 1e6, 1176.45e6
gps_lambda1, gps_lambda2, gps_lambda5 = csts.c/gps_f1, csts.c/gps_f2, csts.c/gps_f5
gps_alpha = gps_f1**2*gps_f2**2/(gps_f1**2-gps_f2**2)/40.318

glonass_f1,glonass_f2 = 1602 * 1e6, 1246 * 1e6

channel_glonass = {"R01":1,"R02":-4,"R03":5,"R04":6,"R05":1,"R06":-4,"R07":5,"R08":6,
           "R09":-2,"R10":-7,"R11":0,"R12":-1,"R13":-2,"R14":-7,"R15":0,"R16":-1,
           "R17":4,"R18":-3,"R19":3,"R20":2,"R21":4,"R22":-3,"R23":3,"R24":2,
            "R25":-5, "R26":-6        
          }

def getFrequency(sv,chanel):
    if sv[0]=='G':
        if chanel[1]=='1': return 1575.42 * 1e6
        elif chanel[1]=='2': return 1227.60 * 1e6
        elif chanel[1]=='5': return 1176.45 * 1e6
    elif sv[0]=='R':
        if chanel[1]=='1': return glonass_f1 + 1e6*channel_glonass[sv]*9/16
        elif chanel[1]=='2': return glonass_f2 + 1e6*channel_glonass[sv]*7/16
    elif sv[0]=='C':
        if chanel[1]=='2': return 1561.098 * 1e6
        elif chanel[1]=='7': return 1207.14 * 1e6
        elif chanel[1]=='6': return 1268.52 * 1e6
    elif sv[0]=='E':
        if chanel[1]=='1': return 1575.42 * 1e6
        elif chanel[1]=='2': return 1227.60 * 1e6
        elif chanel[1]=='5': return 1176.45 * 1e6
        elif chanel[1]=='6': return 1278.75 * 1e6
        elif chanel[1]=='7': return 1207.140 * 1e6
        elif chanel[1]=='8': return 1191.795 * 1e6
    elif sv[0]=='J':
        if chanel[1]=='1': return 1575.42 * 1e6
        if chanel[1]=='2': return 1227.60 * 1e6
        elif chanel[1]=='5': return 1176.45 * 1e6
    elif sv[0]=='S':
        if chanel[1]=='1': return 1575.42 * 1e6
        elif chanel[1]=='5': return 1176.45 * 1e6
    return float('NaN')


def getAlpha(sv,chanel1,chanel2):
    f1 = getFrequency(sv,chanel1)
    f2 = getFrequency(sv,chanel2)
    if f1==float('NaN'): return float('NaN')
    if f2==float('NaN'): return float('NaN')
    if f1==f2: return float('NaN')
    return f1**2*f2**2/(f1**2-f2**2)/40.318
