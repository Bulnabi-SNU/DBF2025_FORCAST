#Root Chord, Tip Chord  (+, -, 몇씩 증가?)
#AR = 5.45, Taper Ratio = 0.65
main_wing_rootC = []
main_wing_tipC = []
main_wing_span = []
main_wing_area = 0.59

vertic_wing_rootC = []
vertic_wing_tipC = []
vertic_wing_span = []
vertic_wing_area = 0.4

horizon_wing_rootC = []
horizon_wing_tipC = []
horizon_wing_span = []
horizon_wing_area = 0.4

main_wing_arr = []
vertic_wing_arr = []
horizon_wing_arr = []


############################## List 저장 ##############################

def rootC(name, operation, num, val, count):
    if name == "main":
        if operation == "+":
            for item in range(count):
                num = num + val * (item + 1)
                main_wing_rootC.append(num)
        elif operation == "-":
            for item in range(len(main_wing_arr)):
                num = num - val * (item + 1)
                main_wing_rootC.append(num)

    elif name == "vertical":
        if operation == "+":
            for item in range(count):
                num = num + val * (item + 1)
                vertic_wing_rootC.append(num)
        elif operation == "-":
            for item in range(count):
                num = num - val * (item + 1)
                vertic_wing_rootC.append(num)

    elif name == "horizontal":
        if operation == "+":
            for item in range(count):
                num = num + val * (item + 1)
                horizon_wing_rootC.append(num)
        elif operation == "-":
            for item in range(count):
                num = num - val * (item + 1)
                horizon_wing_rootC.append(num)



def tipC(name, operation, num, val, count):
    if name == "main":
        if operation == "+":
            for item in range(count):
                num = num + val * (item + 1)
                main_wing_tipC.append(num)
        elif operation == "-":
            for item in range(count):
                num = num - val * (item + 1)
                main_wing_tipC.append(num)

    elif name == "vertical":
        if operation == "+":
            for item in range(count):
                num = num + val * (item + 1)
                vertic_wing_tipC.append(num)
        elif operation == "-":
            for item in range(count):
                num = num - val * (item + 1)
                vertic_wing_tipC.append(num)

    elif name == "horizontal":
        if operation == "+":
            for item in range(count):
                num = num + val * (item + 1)
                horizon_wing_tipC.append(num)          
        elif operation == "-":
            for item in range(count):
                num = num - val * (item + 1)
                horizon_wing_tipC.append(num)


def span(name, operation, num, val, count):
    if name == "main":
        if operation == "+":
            for item in range(count):
                num = num + val * (item + 1)
                main_wing_span.append(num)            
        elif operation == "-":
            for item in range(count):
                num = num - val * (item + 1)
                main_wing_span.append(num)

    elif name == "vertical":
        if operation == "+":
            for item in range(count):
                num = num + val * (item + 1)
                vertic_wing_span.append(num)
        elif operation == "-":
            for item in range(count):
                num = num - val * (item + 1)
                vertic_wing_span.append(num)

    elif name == "horizontal":
        if operation == "+":
            for item in range(count):
                num = num + val * (item + 1)
                horizon_wing_span.append(num)
        elif operation == "-":
            for item in range(count):
                num = num - val * (item + 1)
                horizon_wing_span.append(num)
#######################################################################


############################# 제약 조건 함수 #############################

#if문 조건
def VaildWing(constraints, wing_area, span, taper_ratio):
    return (constraints <= wing_area 
            and span <= 1.8 
            and 0.62 < taper_ratio < 0.67)

#주날개
def MainWingParam(rootC, tipC, span):
    print("[Main Wing]")   

    if isinstance(rootC, list):
        print("-RootChord")
        for i in range(len(rootC)):
            taper_ratio = tipC / rootC[i]
            constraints = ( rootC[i] + tipC ) / 2 * span / 2
            if VaildWing(constraints, main_wing_area, span, taper_ratio):

                main_wing_arr.append([rootC[i], tipC, span])
            else:
                num = i+1
                print(f"{num}번 제한 조건 충족 X")

        
        print(f"\n최종 리스트 = {main_wing_arr}")

    elif isinstance(tipC, list):
        print("-TipChord")
        for i in range(len(tipC)):
            taper_ratio = tipC[i] / rootC
            constraints = ( rootC + tipC[i] ) / 2 * span / 2
            if VaildWing(constraints, main_wing_area, span, taper_ratio):
                main_wing_arr.append([rootC, tipC[i], span])
            else:
                num = i+1
                print(f"{num}번 제한 조건 충족 X")


        print(f"\n최종 리스트 = {main_wing_arr}")

    elif isinstance(span, list):
        print("-Span")
        for i in range(len(span)):
            taper_ratio = tipC / rootC
            constraints = ( rootC + tipC ) / 2 * span[i] / 2
            if VaildWing(constraints, main_wing_area, span, taper_ratio):
                main_wing_arr.append([rootC, tipC, span[i]])
            else:
                num = i+1
                print(f"{num}번 제한 조건 충족 X")

        print(f"\n최종 리스트 = {main_wing_arr}")

    print("\n")

#수직꼬리날개
def VerticalWingParam(rootC, tipC, span):
    print("[Vertical Wing]")   

    if isinstance(rootC, list):
        print("-RootChord")
        for i in range(len(rootC)):
            taper_ratio = tipC / rootC[i]
            constraints = ( rootC[i] + tipC ) / 2 * span / 2
            if VaildWing(constraints, vertic_wing_area, span, taper_ratio):

                vertic_wing_arr.append([rootC[i], tipC, span])
            else:
                num = i+1
                print(f"{num}번 제한 조건 충족 X")

        
        print(f"\n최종 리스트 = {vertic_wing_arr}")

    elif isinstance(tipC, list):
        print("-TipChord")
        for i in range(len(tipC)):
            taper_ratio = tipC[i] / rootC
            constraints = ( rootC + tipC[i] ) / 2 * span / 2
            if VaildWing(constraints, vertic_wing_area, span, taper_ratio):
                vertic_wing_arr.append([rootC, tipC[i], span])
            else:
                num = i+1
                print(f"{num}번 제한 조건 충족 X")


        print(f"\n최종 리스트 = {vertic_wing_arr}")

    elif isinstance(span, list):
        print("-Span")
        for i in range(len(span)):
            taper_ratio = tipC / rootC
            constraints = ( rootC + tipC ) / 2 * span[i] / 2
            if VaildWing(constraints, vertic_wing_area, span, taper_ratio):
                vertic_wing_arr.append([rootC, tipC, span[i]])
            else:
                num = i+1
                print(f"{num}번 제한 조건 충족 X")

        print(f"\n최종 리스트 = {vertic_wing_arr}")

    print("\n")

#수평꼬리날개
def HorizontalWingParam(rootC, tipC, span):
    print("[Horizontal Wing]")   

    if isinstance(rootC, list):
        print("-RootChord")
        for i in range(len(rootC)):
            taper_ratio = tipC / rootC[i]
            constraints = ( rootC[i] + tipC ) / 2 * span / 2
            if VaildWing(constraints, horizon_wing_area, span, taper_ratio):

                horizon_wing_arr.append([rootC[i], tipC, span])
            else:
                num = i+1
                print(f"{num}번 제한 조건 충족 X")

        
        print(f"\n최종 리스트 = {horizon_wing_arr}")

    elif isinstance(tipC, list):
        print("-TipChord")
        for i in range(len(tipC)):
            taper_ratio = tipC[i] / rootC
            constraints = ( rootC + tipC[i] ) / 2 * span / 2
            if VaildWing(constraints, horizon_wing_area, span, taper_ratio):
                horizon_wing_arr.append([rootC, tipC[i], span])
            else:
                num = i+1
                print(f"{num}번 제한 조건 충족 X")


        print(f"\n최종 리스트 = {horizon_wing_arr}")

    elif isinstance(span, list):
        print("-Span")
        for i in range(len(span)):
            taper_ratio = tipC / rootC
            constraints = ( rootC + tipC ) / 2 * span[i] / 2
            if VaildWing(constraints, horizon_wing_area, span, taper_ratio):
                horizon_wing_arr.append([rootC, tipC, span[i]])
            else:
                num = i+1
                print(f"{num}번 제한 조건 충족 X")

        print(f"\n최종 리스트 = {horizon_wing_arr}")

    print("\n")

#######################################################################
############################## Main Code ##############################

MainWingParam()