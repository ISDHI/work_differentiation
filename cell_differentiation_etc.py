# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.6
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% jupyter={"outputs_hidden": true}
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import itertools
import time
# !pip install jupytext --upgrade
np.random.seed(1234)


# %%
def cell_differentiation(cell_list, limit, probR, probS, probZ,print_option=0, rule1="B", rule2="C", lowcycle=20, highcycle=40):
    
    '''
    細胞を分化、分裂させる関数。
    cell_listには細胞の初期状態リスト、
    rには分裂、分化させる回数、
    print_optionはプリントの詳細度合いの引数。
    probは細胞が分化する確率
    アウトプットの辞書の形式
    {roop数:{"S":Sの数, "A":Aの数, "B":Bの数, "G":Gの数}}のフォーマットの辞書
    '''
    cell_dict = {}
    #print("細胞の初期状態は",cell_list, "です。", "細胞の分化の確率は",probR,"です。", f"細胞の保持の確率は{probS}です",f"細胞が消失する確率は{probZ}です")
    n = 0 #作業の試行回数
    r = -1 #forループの回数
    while r <= limit:
        r += 1
        if r >= limit + 1 :
                    print("最終的な細胞は以下になります", "A",(cell_list.count("A") + (2 * cell_list.count("AA")), rule1,cell_list.count(rule1), rule2, (2 * cell_list.count("CC") + cell_list.count(rule2))))
                    return cell_dict
                    break
        #print(cell_list)
        #print(cell_dict)
        cell_dict[r] = {"S": cell_list.count("S"), "A": (cell_list.count("A") + (2*cell_list.count("AA"))), rule1: cell_list.count(rule1), rule2: (cell_list.count(rule2) + (2*cell_list.count("CC"))), "AA":cell_list.count("AA")}
        joined_cell_str = "".join(cell_list)
        cell_list = list(joined_cell_str)
        if r != 0:
            for c in range(cell_list.count("Z")):
                cell_list.remove("Z")
            
            for i in range(int(cell_list.count(rule1) + np.floor(0.5 * cell_list.count(rule2)))):
                        if "A" in cell_list:
                            cell_list.remove("A")
                        else:
                            pass
            if "B" in cell_list and "C" in cell_list:
                if cell_list.count(rule2) > cell_list.count(rule1):
                    for a in range(cell_list.count(rule2)):
                        if "A" in cell_list:
                               cell_list.remove("A")
                        else:
                               pass
                elif cell_list.count(rule1) > cell_list.count(rule2):
                    for a in range(cell_list.count(rule1)):
                        if "A" in cell_list:
                               cell_list.remove("A")
                        else:
                               pass
                           
                elif cell_list.count(rule2) == cell_list.count(rule1):
                    for a in range(cell_list.count(rule2)):
                        if "A" in cell_list:
                               cell_list.remove("A")
                        else:
                               pass
                               
        trans_list = {"A": ["AA", "A", "B", "CC", "Z"], "S":"A", "B":"Z", "C":"Z"}
        
        if r >= limit + 1 :
                    #print("最終的な細胞は以下になります", "A",cell_list.count("A"), rule1,cell_list.count(rule1), rule2, cell_list.count(rule2))
                    cell_dict[r] = {"S":cell_list.count("S"), "A":(cell_list.count("A") + cell_list.count("AA")), rule1:cell_list.count(rule1), rule2:(cell_list.count(rule2) + cell_list.count("CC")), "AA":cell_list.count("AA")}
                    return cell_dict
                    break
        if len(cell_list) == 0:
            print("細胞は消失しました", "ループ数", r)
            print(cell_dict)
            return cell_dict
            r += 10000
            break
        elif r > 0 and r <= lowcycle and cell_list.count("A") >= 10000:
            #print(f"細胞サイクルが{lowcycle}にもかかわらず細胞Aの数が10000を超えたので停止しました。")
            #print(r)
            return cell_dict
            break
        elif r > 0 and cell_list.count("A") >= 15000:
            #print(cell_list.count("A"))
            #print("細胞Aが15000を超えたので停止しました。")
            #print(r)
            return cell_dict
            break
        elif r > 0 and r >= highcycle and cell_list.count("A") <= 3 and cell_list.count("C") <= 3 and cell_list.count("B") <= 3 :
            #print(f"細胞サイクルが{highcycle}にもかかわらず、細胞Aの数が<=20のため停止しました。")
            #print(r)
            return cell_dict
            break
        elif r == 0:
            #print("forループは終了しました。新たな細胞分裂、分化を始めます。", cell_list, "この細胞群を使用します", "現在のループ回数",r, "現在の細胞数",len(cell_list))
            for index, cell in enumerate(cell_list):
                #index = np.random.choice(range(len(cell_list_test)))
                cell_type = cell
                if cell_type == "A":
                    cell_type_for_difference = trans_list[cell_type]
                    probabilities = [1 - (2*probR) - probS - probZ, probS, probR, probR, probZ]
                    chosen_one = np.random.choice(cell_type_for_difference, p=probabilities)
                    cell_list[index] = chosen_one
                    if chosen_one == "Z":
                        cell_list.remove("Z")
                    else:
                        pass
                    n += 1
                    if print_option == 1:
                        print("作業後の細胞の状態",cell_list, "細胞を分化、分裂させました", "作業対象の細胞",cell, index + 1)
                elif cell_type == "B" or cell_type == "C":
                    cell_list[index] = trans_list[cell_type]
                    n += 1
                    if print_option == 1:
                        print("作業後の細胞の状態",cell_list, "細胞を老化させました", "作業対象の細胞",cell, index + 1)
                #elif cell_type == "B" or cell_type == "C":
                    #cell_list.remove(cell_type)
                    #n += 1
                    #if print_option == 1:
                        #print("作業後の細胞の状態",cell_list, "老化した細胞を取り除きました", "作業対象の細胞",cell, index + 1)
                elif cell_type == "S":
                    cell_list.append("A")
                    n+= 1
                    if print_option == 1:
                        print("作業後の細胞の状態",cell_list)
        else:
            #for i in range(int(cell_list.count(rule1) + cell_list.count(rule2))):
                #if "A" in cell_list:
                    #cell_list.remove("A")
                #else:
                    #pass
                    
            for index, cell in enumerate(cell_list):
                #index = np.random.choice(range(len(cell_list_test)))
                cell_type = cell
                if cell_type == "A":
                    cell_type_for_difference = trans_list[cell_type]
                    probabilities = [1 - (2*probR) - probS - probZ, probS, probR, probR, probZ]
                    chosen_one = np.random.choice(cell_type_for_difference, p=probabilities)
                    cell_list[index] = chosen_one
                    if chosen_one == "Z":
                        cell_list.remove("Z")
                    else:
                        pass                    
                    n += 1
                    if print_option == 1:
                        print("作業後の細胞の状態",cell_list, "細胞を分化、分裂させました", "作業対象の細胞",cell, index + 1)
                elif cell_type == "B" or cell_type == "C":
                    cell_list[index] = trans_list[cell_type]
                    n += 1
                    if print_option == 1:
                        print("作業後の細胞の状態",cell_list, "細胞を老化させました", "作業対象の細胞",cell, index + 1)
                #elif cell_type == "B" or cell_type == "C":
                    #cell_list.remove(cell_type)
                    #n += 1
                    #if print_option == 1:
                        #print("作業後の細胞の状態",cell_list, "老化した細胞を取り除きました", "作業対象の細胞",cell, index + 1)
                elif cell_type == "S":
                    cell_list.append("A")
                    n+= 1
                    if print_option == 1:
                        print("作業後の細胞の状態",cell_list)


# %%
def cell_dict_toplot(cell_dict):
    #print("time, cell_A_list, cell_B_list, cell_G_list　の順番でリストを１つのタプルとして返します")
    time = []
    cell_A_list = []
    cell_B_list = []
    cell_C_list = []
    '''
    インプットは{roop数:{"S":Sの数, "A":Aの数, "B":Bの数, "G":Gの数}}のフォーマットの辞書
    '''
    
    for keys, values in cell_dict.items():
        time.append(keys)
        for key, value in values.items():
            if key == "A":
                cell_A_list.append(value)
            elif key == "B":
                cell_B_list.append(value)
            elif key == "C":
                cell_C_list.append(value)
    return time, cell_A_list, cell_B_list, cell_C_list


# %%
def plot_cellfig(cell_lists, sprob, rprob, tprob):
    '''
    引数はcell_dict_toplotによって作成された
    time, cell_A_list, cell_B_list, cell_G_list　のタプル
    細胞A,B,Gの時系列変化と細胞B,Gの相図を作成する関数。
    '''
    time = cell_lists[0]
    cell_list_A = cell_lists[1]
    cell_list_B = cell_lists[2]
    cell_list_C = cell_lists[3]
    
    if len(cell_list_A) >= 150:#何サイクル以上回すことができたら図をプロットするのか記述
        print(f"確率はs:{sprob}、r:{rprob}、t:{tprob}です")
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(15, 7 * len(cell_lists)))
        ax1.plot(time, cell_list_A, color="green")
        ax1.set_title('A_count')
        ax1.set_xlabel('time')
        ax1.set_ylabel('cell_count_A')
    
        ax2.plot(time, cell_list_B, color="green")
        ax2.set_title('B_count')
        ax2.set_xlabel('time')
        ax2.set_ylabel('cell_count_B')
    
        ax3.plot(time, cell_list_C, color="green")
        ax3.set_title('C_count')
        ax3.set_xlabel('time')
        ax3.set_ylabel('cell_count_C')
    
        ax4.plot(cell_list_A, cell_list_C, color="green")
        ax4.set_title('plot between A and C')
        ax4.set_xlabel("cell_count_A")
        ax4.set_ylabel("cell_count_C")
        plt.tight_layout()
        plt.show()
    else:
        pass


# %%
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def phase_diagram(suc_S, suc_R, miss_S, miss_R, suc_T, miss_T):
    # sとrの範囲を示すグラフ
    x = np.linspace(0, 1, 100)
    y = 1 - (2 * x)

    # 3次元プロットを作成
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # x, y, z軸のラベルとタイトルを設定
    ax.set_xlabel('Rptob')
    ax.set_ylabel('Sprob')
    ax.set_zlabel('T_prob')
    ax.set_title('Phase Diagram')

    # 3次元散布図を作成
    ax.scatter(suc_R, suc_S, suc_T, color='blue', s=5)
    ax.scatter(miss_R, miss_S, miss_T, color='red', s=5)

    # x, y, z軸の範囲を設定
    ax.set_xlim(0, 0.6)
    ax.set_ylim(0, 1)

    # グリッドを表示
    ax.grid(True)

    # プロットを表示
    plt.show()
    


# %%
#お試し用のセル
np.random.seed(12345)
cell_dict = cell_differentiation(["S"], 2000, 0.07, 0.02, 0.2,print_option=0, )
#print(cell_dict)
#print(len(cell_dict))
cell_lists = cell_dict_toplot(cell_dict)
plot_cellfig(cell_lists, 0.07, 0.11, 0)

# %%
#probZ = 0.1用

results = []

# sとrの範囲を設定

s_values = np.arange(0, 1.01, 0.01) 
r_values = np.arange(0, 1.01, 0.01)


#for s, r in zip(s_values, rvalues:
    #results.append((s, r))
# グリッドサーチ
for s, r in itertools.product(s_values, r_values):
    s = round(s, 2)
    r = round(r, 2)
    if  (1 - (2*r) - s) + s + r + r == 1 and s != 0 and r != 0 and (1 - (2*r) - s) != 0 and r*2 <= 1 and s + (r*2) <=1 and (1 - (2*r) - s) > 0 and (1 - (2*r) - s) <= 1 and (1 - (2*r) - s) > 0.1:
        results.append((s, r))

# 結果の表示
#for s, r in results:
    #print(f"s: {s}, r: {r}")
    
    
#del results[0]
#print(results)

# %%
'''
suc_S = []
suc_R = []
miss_S = []
miss_R = []

start_time = time.time()

for s, r in results:
    np.random.seed(1234)
    cell_dict = cell_differentiation(["S"], 200, r, s,0)
    
    if len(cell_dict) >= 150:
        suc_S.append(s)
        suc_R.append(r)
    else:
        miss_S.append(s)
        miss_R.append(r)
    cell_lists = cell_dict_toplot(cell_dict)
    plot_cellfig(cell_lists, s, r)
    
phase_diagram(suc_S, suc_R, miss_S, miss_R)

end_time = time.time()
process_time = end_time - start_time
print(process_time)
'''


# %%
#phase_diagram(suc_S, suc_R, miss_S, miss_R)

# %%
prob_Z = [i / 10 for i in range(0, 10, 1)]
print(prob_Z)

# %%
"""
for s, r in results:
    np.random.seed(1234)
    for z in prob_Z:
        if np.float(s + 2*r + z) <= 0.9:
            cell_dict = cell_differentiation(["S"], 200, r, s,z)
            cell_lists = cell_dict_toplot(cell_dict)
            plot_cellfig(cell_lists)
            
"""

# %%
t_values = np.arange(0, 0.31, 0.05) 

s_values = np.arange(0, 1.01, 0.01) 
r_values = np.arange(0, 1.01, 0.01)

print(len(s_values))

param = []


for t in t_values:
    m = 0
    n = 0
    for s in s_values:
        for r in r_values:
            if  (1 - (2*r) - s - t) + s + r + r + t == 1 and s != 0 and r != 0 and t!= 0 and (1 - (2*r) - s - t) != 0 and r*2 <= 1 and s + (r*2) + t <=1 and (1 - (2*r) - s - t) > 0 and (1 - (2*r) - s - t) <= 1 and (1 - (2*r) - s - t) > 0.1:
                    param.append((s, r, t))
print(len(param))

# %% jupyter={"outputs_hidden": true}
suc_S = []
suc_R = []
miss_S = []
miss_R = []
suc_T = []
miss_T = []
start_time = time.time()

for s, r, t in param:
    np.random.seed(1234)
    cell_dict = cell_differentiation(["S"], 200, r, s, t)#r,s,tの順番注意！
    
    if len(cell_dict) >= 150:
        suc_S.append(s)
        suc_R.append(r)
        suc_T.append(t)
    else:
        miss_S.append(s)
        miss_R.append(r)
        miss_T.append(t)
    cell_lists = cell_dict_toplot(cell_dict)
    plot_cellfig(cell_lists, s, r, t)
    
phase_diagram(suc_S, suc_R, miss_S, miss_R, suc_T, miss_T)

end_time = time.time()
process_time = end_time - start_time
print(process_time)

# %%
