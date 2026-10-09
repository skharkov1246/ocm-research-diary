import importlib.util,sys
spec=importlib.util.spec_from_file_location('c','check.py')
import io,contextlib
with contextlib.redirect_stdout(io.StringIO()):
    c=importlib.util.module_from_spec(spec); spec.loader.exec_module(c)
P=c.project
def sc_run(sc,capex0=12.0,staff0=0.6,olig=0.85,price=3200,decl=0.0,real=True,gasp=0.25):
    p=dict(Y=0.35,olig=olig,price=price,capex=capex0*sc**0.7,gas=0,var=250,staff=staff0*sc**0.4,sec=0.3*sc**0.3,logi=150,maint=0.035,
           ch4=(c.ch4_t if real else 7300)*sc,decl=decl)
    if real: p.update(avail=0.9,ramp=0.6,tax=0.34,wc=0.2,owner=0.1)
    r=P(**p)
    # subtract gas cost properly scaled (approx: undiscounted annuity adj via e1 only reported)
    return r
for sc in (5,15):
    print(sc,'their A        ',sc_run(sc,real=False))
    print(sc,'+mech realism  ',sc_run(sc))
    print(sc,'+olig .75      ',sc_run(sc,olig=0.75))
    print(sc,'+staff 1.2     ',sc_run(sc,olig=0.75,staff0=1.2))
    print(sc,'+price 2600    ',sc_run(sc,olig=0.75,staff0=1.2,price=2600))
    print(sc,'+capex0 22     ',sc_run(sc,olig=0.75,staff0=1.2,price=2600,capex0=22))
    print(sc,'+decl 10%      ',sc_run(sc,olig=0.75,staff0=1.2,price=2600,capex0=22,decl=0.10))
    print(sc,'decene $3200,capex22,olig.75,real,no decl',sc_run(sc,olig=0.75,staff0=1.2,capex0=22))
