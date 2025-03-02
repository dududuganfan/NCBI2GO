from Bio import SeqIO
import pandas as pd
import re
from taxonomy_ranks import TaxonomyRanks

# 读取GenBank文件, 设置Excel文件存储位置
gb_file = "DataCleaning/Aves_RefSeq.gb"
excel_file = "DataCleaning/Aves_RefSeq.xlsx"
records = SeqIO.parse(gb_file, "genbank")

# 读取规范化字典(如果没有字典，则设置一个空字典)
dict_df = pd.read_excel("DataCleaning/Dictionary.xlsx")

# 创建一个空的数据帧,对应LOCUS、带坐标的Product、不带坐标调零的Product等等
data = {
    'LOCUS': [],
    'Superorder': [],
    'Order': [],
    'Suborder': [],
    'Superfamily': [],
    'Family': [],
    'Subfamily': [],
    'Genus': [],
    'Subgenus': [],
    'Species_name': [],
    'circular/linear': [],
    'completeness': [],
    '基因总长': [],
    'Product': [],
    'Product_No_Location': [],
    }


# 针对字典的规范化的函数,从字典中查询是否与value值相等，若相等则替换为new_value,否则不变
def normalize_value(normal_value):
    replacement = dict_df[dict_df['value'].str.lower() == normal_value.lower()]['new_value'].values
    if len(replacement) > 0:
        return replacement[0]
    else:
        return normal_value


# 处理Product字符串反向的情况,进行翻转并调整符号
def process_product_Fend(product_string):
    # 以@拆分字符串元素
    components = product_string.split('@')
    # 存储字符串元素的数组
    processed_components = []
    # 初始化结果
    processed_product: str = ""
    # 从最后一个元素开始遍历,反向存入数组
    for i in range(len(components) - 1, -1, -1):
        # 如果该元素存在-则删除
        if components[i]:
            if components[i][0] == '-':
                components[i] = components[i][1:]
            # 如果不存在-则添加
            else:
                components[i] = '-' + components[i][0:]
            # 添加到处理后的数组
        processed_components.append(components[i])
        # 重新连接为字符串，以@进行分隔
        processed_product = '@'.join(processed_components)
    return processed_product


# 遍历GenBank记录并提取位置信息及相应的值
for record in records:
    # 提取LOCUS
    locus = record.name
    # 提取线形or环形
    circular = record.annotations.get("topology", "linear")
    # 获取GenBank文件中的定义行（DEFINITION）
    definition_line = record.description
    completeness = ""
    if "partial" in definition_line:
        completeness = "partial"
    elif "complete" in definition_line:
        completeness = "complete"
    # 如数据为NC版本,需要提取来源gene_source,初始化为空, gene_source存储NC数据来源的LOCUS,最后需要将NC与来源合并,使用/分隔
    gene_source = ""
    if locus.startswith("NC_"):
        comment = record.annotations["comment"]
        # 来源一般跟随在 from 或者 identical to 之后,提取这两个位置之后的字符串获取来源LOCUS
        if "from " in comment:
            gene_source = comment.split("from ")[1].split(".")[0]
        if "identical to " in comment:
            gene_source = comment.split("identical to ")[1].split(".")[0]
        # 将NC和来源LOCUS合并,中间使用/进行分隔
        locus = gene_source + "/" + locus

    # 提取物种名及阶元(使用taxonomy_ranks中的TaxonomyRanks,其来源为ETE3,该库会从NCBI的ftp服务器下载整个分类数据库)
    species_name = record.annotations["organism"]
    # rank_taxon = TaxonomyRanks(species_name)
    # if rank_taxon:
    #     rank_taxon.get_lineage_taxids_and_taxanames()
    # # 分别提取八个阶元,如果为空值,会自动填入NA
    # ranks = ('superorder', 'order', 'suborder', 'superfamily', 'family', 'subfamily', 'genus', 'subgenus')
    # # 将提取出来的阶元存入数组中
    # result_array = []
    # for potential_taxid in rank_taxon.lineages:
    #     for rank in ranks:
    #         taxon, taxid_of_taxon = rank_taxon.lineages[potential_taxid][rank]
    #         result_array.append(taxon)
    try:
        rank_taxon = TaxonomyRanks(species_name)
        if rank_taxon:
            rank_taxon.get_lineage_taxids_and_taxanames()
            # Define the ranks to be extracted
            ranks = ('superorder', 'order', 'suborder', 'superfamily', 'family', 'subfamily', 'genus', 'subgenus')
            # Initialize the result array
            result_array = []
            for potential_taxid in rank_taxon.lineages:
                for rank in ranks:
                    taxon, taxid_of_taxon = rank_taxon.lineages[potential_taxid].get(rank, ("NA", "NA"))
                    result_array.append(taxon)
    except ValueError as e:
        print(f"Skipping species {species_name} due to error: {e}")

    # 依次将数组中的值赋值给各个阶元变量
    subgenus = result_array[7]
    genus = result_array[6]
    subfamily = result_array[5]
    family = result_array[4]
    superfamily = result_array[3]
    suborder = result_array[2]
    order = result_array[1]
    superorder = result_array[0]

    # 创建一个products的集合[location,product]
    products = {}
    # 设置一个last_start_position 和 last_end_position 用于判断包含关系
    last_value = None
    last_end_location = None
    last_start_position = None

    # 初始化feature及gene_length
    feature: object
    gene_length: str = ""
    # 遍历该记录中的所有location
    for feature in record.features:
        # 排除基因的总长位置,并通过该Location提取基因长度.如果type = source 则跳出循环，进行下一次循环
        if feature.type == "source":
            start_position = feature.location.start + 1
            end_position = feature.location.end
            # gene_length 只取末尾
            # gene_length = f"{start_position}..{end_position}"
            gene_length = end_position
            continue

        # 依次遍历注释信息
        # 1、定位location
        # 获取location的start_position 和 end_position
        start_position = feature.location.start + 1
        end_position = feature.location.end
        # 获取strand(+省略)
        strand = "" if feature.strand == 1 else "-"
        # 使用正则去除start_position和end_position中的非数字字符(非数字如:><符号)
        start_position = re.sub(r'\D', '', str(start_position))
        end_position = re.sub(r'\D', '', str(end_position))
        # 对location进行赋值,格式:start_position..end_position(strand)
        location = f"{start_position}..{end_position}({strand})"
        # 如果存在start = end 的情况(即单个基因位置),则直接跳过该循环
        if start_position == end_position:
            continue

        # value 初始化,用于存储每个location对应的具体注释信息,通过系列判断对其进行赋值
        value: str = ""
        # 2、提取type进行分析，对不同type进行不同处理
        # (2.1) type = CR 的情况:
        if feature.type in ["D-loop", "repeat_region", "rep_origin", "regulatory",
                            "C_region", "stem_loop", "N_region", "misc_RNA", "Region"]:
            value = "CR"
        # (2.2) type = GAP 的情况:
        elif feature.type.lower() in ["gap"] or "gap" in feature.type.lower():
            value = "GAP"
        # (2.3) type = tRNA\rRNA\CDS的情况
        elif feature.type in ["tRNA", "rRNA", "CDS"]:
            # 优先级顺序: product>gene>note>else
            # rRNA
            if feature.type == "rRNA":
                # 提取 product or gene 中的值
                if "product" in feature.qualifiers:
                    value = feature.qualifiers["product"][0]
                elif "gene" in feature.qualifiers:
                    value = feature.qualifiers["gene"][0]
                elif "note" in feature.qualifiers:
                    value = feature.qualifiers["note"][0]
                # rrnS 的情况
                if "12" in value.lower() or "rrns" in value.lower() or "rns" in value.lower() or value.lower().startswith('s'):
                    value = "rrnS"
                # rrnL 的情况
                elif "16" in value.lower() or "rrnl" in value.lower() or "rnl" in value.lower() or value.lower().startswith('l'):
                    value = "rrnL"
                else:
                    value = "?"
            # tRNA
            elif feature.type == "tRNA":
                # 提取 product or gene 中的值
                if "product" in feature.qualifiers:
                    value = feature.qualifiers["product"][0]
                elif "gene" in feature.qualifiers:
                    value = feature.qualifiers["gene"][0]
                elif "note" in feature.qualifiers:
                    value = feature.qualifiers["note"][0]
                # A
                if "ala" in value.lower() or "trnA" in value:
                    value = "A"
                # C
                elif "cys" in value.lower() or "trnC" in value:
                    value = "C"
                # D
                elif "asp" in value.lower() or "trnD" in value:
                    value = "D"
                # E
                elif "glu" in value.lower() or "trnE" in value:
                    value = "E"
                # F
                elif "phe" in value.lower() or "trnF" in value:
                    value = "F"
                # G
                elif "gly" in value.lower() or "trnG" in value:
                    value = "G"
                # H
                elif "his" in value.lower() or "trnH" in value:
                    value = "H"
                # I
                elif "ile" in value.lower() or "trnI" in value:
                    value = "I"
                # K
                elif "lys" in value.lower() or "trnK" in value:
                    value = "K"
                # L(区分L1和L2)
                elif "leu" in value.lower() or "trnl" in value.lower() or "trna-l" in value.lower():
                    value = "L"
                # M
                elif "met" in value.lower() or "trnM" in value:
                    value = "M"
                # N
                elif "asn" in value.lower() or "trnN" in value:
                    value = "N"
                # P
                elif "pro" in value.lower() or "trnP" in value:
                    value = "P"
                # Q
                elif "gln" in value.lower() or "trnQ" in value:
                    value = "Q"
                # R
                elif "arg" in value.lower() or "trnR" in value:
                    value = "R"
                # S(区分S1和S2)
                elif "ser" in value.lower() or "trns" in value.lower() or "trna-s" in value.lower():
                    value = "S"
                # T
                elif "thr" in value.lower() or "trnT" in value:
                    value = "T"
                # V
                elif "val" in value.lower() or "trnV" in value:
                    value = "V"
                # W
                elif "trp" in value.lower() or "trnW" in value:
                    value = "W"
                # Y
                elif "tyr" in value.lower() or "trnY" in value:
                    value = "Y"
                else:
                    value = "?"
            # CDS
            elif feature.type == "CDS":
                # 提取 product or gene 中的值
                if "product" in feature.qualifiers:
                    value = feature.qualifiers["product"][0]
                elif "gene" in feature.qualifiers:
                    value = feature.qualifiers["gene"][0]
                elif "note" in feature.qualifiers:
                    value = feature.qualifiers["note"][0]
                # cob
                if "cytochrome b" in value.lower() or "cytochrome-b" in value.lower() \
                        or "cytochtrome b" in value.lower() or "cytb" in value.lower() \
                        or "cyt b" in value.lower() or "ctyb" in value.lower():
                    value = "cob"
                # cox
                elif "cytochrome c" in value.lower() or "cytochrome o" in value.lower() or "CO" in value:
                    if "3" in value or "III" in value:
                        value = "cox3"
                    elif "2" in value or "II" in value:
                        value = "cox2"
                    elif "1" in value or "I" in value:
                        value = "cox1"
                    else:
                        value = "?"
                # atp
                elif "atp" in value.lower() or "triphosphatase" in value.lower():
                    if "6" in value:
                        value = "atp6"
                    elif "8" in value:
                        value = "atp8"
                    else:
                        value = "?"
                # nad
                elif "nicotinamide" in value.lower() or "nd" in value.lower() or "nadh" in value.lower()\
                        or "nad" in value.lower():
                    if "3" in value or "III" in value:
                        value = "nad3"
                    elif "2" in value or "II" in value:
                        value = "nad2"
                    elif "4" in value or "IV" in value:
                        if "4l" in value.lower() or "4 l" in value.lower() \
                                or "ivl" in value.lower() or "iv l" in value.lower():
                            value = "nad4l"
                        else:
                            value = "nad4"
                    elif "6" in value or "VI" in value:
                        value = "nad6"
                    elif "5" in value or "V" in value:
                        value = "nad5"
                    elif "1" in value or "I" in value:
                        value = "nad1"
                    else:
                        value = "?"
                else:
                    value = "?"
        # (2.4) type = misc_feature/misc_.../gene的情况
        elif feature.type in ["misc_feature", "gene"]:
            if 'product' in feature.qualifiers:
                # 提取product中的值
                value = feature.qualifiers['product'][0]
            elif "gene" in feature.qualifiers:
                # 提取gene中的值
                value = feature.qualifiers['gene'][0]
            elif "note" in feature.qualifiers:
                value = feature.qualifiers['note'][0]
            else:
                # tRNA\rRNA\CDS\gene中不包含product、gene属性(防止特殊情况)
                value = " ".join([f"{feature.type} Key: {key}, Value: {value}"
                                  for key, value in feature.qualifiers.items()])
            value = normalize_value(value)  # 基于字典规范化处理
            # 如果包含region或者repeat,则将value赋值为CR
            if "region" in value or "repeat" in value or "CSB" in value\
                    or "loop" in value or "CR" in value or "replication" in value:
                value = "CR"
        # (2.5) type = mRNA/unsure 的情况 则替换名称为?,进mitos,在鸟类中mRNA/unsure视为错误
        elif feature.type in ["mRNA", "unsure"]:
            value = "?"
        # (2.6) type = else (exon,intron,variation,misc_difference)(2和4处理相同,可以进行合并处理,后续看有无4新的处理方式)
        else:
            if 'product' in feature.qualifiers:
                # 提取product中的值
                value = feature.qualifiers['product'][0]
            elif "gene" in feature.qualifiers:
                # 提取gene中的值
                value = feature.qualifiers['gene'][0]
            elif "note" in feature.qualifiers:
                value = feature.qualifiers['note'][0]
            else:
                # tRNA\rRNA\CDS\gene中不包含product、gene属性(防止特殊情况)
                value = " ".join([f"{feature.type} Key: {key}, Value: {value}"
                                  for key, value in feature.qualifiers.items()])
            value = normalize_value(value)  # 基于字典规范化处理
            # 如果包含region或者repeat,则将value赋值为CR
            if "region" in value or "A-T Rich Region" in value or "repeat" in value or "CSB" in value\
                    or "loop" in value or "CR" in value or "replication" in value:
                value = "CR"

        if "Key" in value:
            continue
        # 将location和value存入products中
        products.setdefault(location, value)

        # 判断同一个location在products中是否已经存在值，如果存在值需要进行判断(使用前值还是后值?)
        if location in products:
            # 当同一个location对应多个value时,优先保留tRNA rRNA CDS 中的注释值
            if value != products[location] and feature.type not in ["tRNA", "rRNA", "CDS"]:
                products[location] = last_value
            else:
                products[location] = value

        # 3、判断当前location是否与上一个location存在包含关系——(判断包含存在逻辑问题待解决!!!)
        if last_end_location is not None and last_start_position is not None and (value == last_value or value == "?"):
            if last_end_location > end_position and last_start_position <= start_position:
                del products[location]
                end_position = last_end_location
                start_position = last_start_position
                value = last_value
            elif last_start_position < start_position and last_end_location >= end_position:
                del products[location]
                end_position = last_end_location
                start_position = last_start_position
                value = last_value

        # 循环结束之前将当前值赋值给last,用于下一次循环
        last_start_position = start_position
        last_end_location = end_position
        last_value = value

    # 将所有注释信息合并为一个字符串，使用@分隔
    product_str = ''.join([f"@{location}@{value}" for location, value in products.items()])
    # str 用来存放包含location的product, str1用来存放不包含location的product
    strand_pattern = r"\((.*?)\)"  # 编写正则从location中提取符号(+/-)
    product_str1 = '@'.join([f"{re.search(strand_pattern, location).group(1)}{value}"
                             for location, value in products.items()])

    # 判断product_str1是否符合规范，如果不是则进行调零or翻转
    index = product_str1.find('F')
    if index == len(product_str1) - 1 and product_str1.endswith("-F"):
        # 情况1:F在最后一个(len(product_str1) - 1) 将整个字符串进行翻转,并改变正负号
        product_str1 = process_product_Fend(product_str1)
    elif index not in [0, -1]:
        # 情况: F不在开头,调零
        if "-F" in product_str1:
            product_str1 = product_str1[index + 2:] + '@' + product_str1[:index + 1]
            product_str1 = process_product_Fend(product_str1)
        else:
            product_str1 = product_str1[index:] + '@' + product_str1[:index - 1]
    else:
        # 保持不变
        product_str1 = product_str1

    # 将所有的LOCUS 和 Product 存入data
    data['LOCUS'].append(locus)
    data['Superorder'].append(superorder)
    data['Order'].append(order)
    data['Suborder'].append(suborder)
    data['Superfamily'].append(superfamily)
    data['Family'].append(family)
    data['Subfamily'].append(subfamily)
    data['Genus'].append(genus)
    data['Subgenus'].append(subgenus)
    data['Species_name'].append(species_name)
    data['circular/linear'].append(circular)
    data['completeness'].append(completeness)
    data['基因总长'].append(gene_length)
    data['Product'].append(product_str)
    data['Product_No_Location'].append(product_str1)

# 处理NC数据与来源数据的合并,并列出两者注释信息之间的不同
# 创建一个数据帧,并将LOCUS进行分列
df = pd.DataFrame(data)

df[['LOCUS_1', 'LOCUS_2']] = df['LOCUS'].str.split('/', n=-1, expand=True)

# 创建一个新列'Remark',初始化为空字符串,比较并添加备注
df['Remark'] = ''
# 比较NC与来源的product不同之处
for index, row in df.iterrows():
    duplicate_rows = df[df['LOCUS_1'] == row['LOCUS_1']]
    if len(duplicate_rows) > 1:
        for col in ['Product_No_Location']:
            unique_values = duplicate_rows[col].str.split('@').dropna().tolist()
            if len(set(map(tuple, unique_values))) > 1:
                # 列举所有不同的地方
                diff_indices = []
                for values in zip(*unique_values):
                    if len(set(values)) > 1:
                        diff_indices.append(f"[{values[1]}:{values[0]}]")
                if diff_indices:
                    df.at[index, 'Remark'] = f"{','.join(diff_indices)}"

# 添加原始索引列(用于合并)
df['Original_Index'] = df.index
# 根据LOCUS_1进行去重，保留LOCUS_2不为空的那一列
df = df.sort_values(by=['LOCUS_1', 'LOCUS_2'], na_position='last').drop_duplicates(subset='LOCUS_1', keep='first')
# 恢复原始顺序
df = df.sort_values(by='Original_Index').drop('Original_Index', axis=1)

# 合并LOCUS_1和LOCUS_2列，如果LOCUS_2中没有值则使用LOCUS_1的值
df['LOCUS'] = df.apply(
    lambda row_data:
    row_data['LOCUS_2'] + '/' + row_data['LOCUS_1']
    if pd.notnull(row_data['LOCUS_2'])
    else row_data['LOCUS_1'], axis=1
)

# 计算Product_No_Location列的计数并添加到新列Count
df['Count'] = df['Product_No_Location'].map(df['Product_No_Location'].value_counts())

# 将数据帧保存为Excel文件
df.to_excel(excel_file, index=False)
