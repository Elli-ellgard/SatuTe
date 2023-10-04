import pandas as pd
import os
from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import glob
import argparse
import re

parser = argparse.ArgumentParser(description="Nodestyle")
parser.add_argument("-image", action="store_true", help="Enable image generation")
parser.add_argument("-midpoint", action="store_true", help="Use midpoint as outgroup")
parser.add_argument("-cycle", action="store_true", help="creates image in cycle form")
parser.add_argument("-compact", action="store_true", help="makes image more compact")

args = parser.parse_args()
directory_path = os.environ.get("SATUTE_DIRECTORY_PATH")

if directory_path is None:
    print("Error: SATUTE_DIRECTORY_PATH environment variable not set.")
    exit(1)

image_directory = os.path.join(directory_path, "image")  # New directory path

if not os.path.exists(image_directory):
    os.makedirs(image_directory)

treefile_path = glob.glob(os.path.join(directory_path, f"*{os.path.extsep}node.treefile"))[0]






pattern = os.path.join(directory_path, f"*{os.path.extsep}treefile")
treefile_paths = glob.glob(pattern)

# Filter out the paths with ".node.treefile"
filtered_treefile_paths = [path for path in treefile_paths if not path.endswith(".node.treefile")]

print("Matching .treefile paths:")
print(filtered_treefile_paths)

files = os.listdir(directory_path)



matching_csv_files = [file for file in files if file.endswith(".csv")]


if filtered_treefile_paths:
    treefile_default_path = filtered_treefile_paths[0]

    with open(treefile_default_path, "r") as treefile:
        treefile_content = treefile.read()

    pattern = r'\)(\d+):'

    matches = re.findall(pattern, treefile_content)

    for csv_file_name in matching_csv_files:
        csv_file_path = os.path.join(directory_path, csv_file_name)
        df2 = pd.read_csv(csv_file_path)
        edge_column = df2['edge']
        df2['Node_A'] = edge_column.str.extract(r'\(([^,]+),', expand=False)
        df2['Node_B'] = edge_column.str.extract(r',([^,)]+)', expand=False)

        with open(treefile_path, 'r') as f:
            s = f.read()

        last_occurrence_index = s.rfind(';')

        s = s[:last_occurrence_index] + 'Node1' + s[last_occurrence_index:]

        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_colwidth', None)

        df2 = df2[
            (df2['Node_B'].str.strip().str.startswith("Node")) &
            (df2['Node_A'].str.strip().str.startswith("Node"))
            ]

        df2['Node_A-Node_B'] = df2.apply(lambda row: (row['Node_A'] + '-' + row['Node_B']).strip(), axis=1)
        #df2['Node_A-Node_B'] = df2['Node_A-Node_B'].str.replace('*', '', regex=False)
        df2['Node_A-Node_B'] = df2['Node_A-Node_B'].str.replace(r'\s+', '', regex=True)

        if matches:
            #print(matches)
            def parentheses_function(s):
                stack = []
                parentheses = []
                for i, char in enumerate(s):
                    if char == '(':
                        stack.append(i)
                    elif char == ')':
                        if stack:
                            open_idx = stack.pop()
                            parentheses.append((open_idx, i))
                return parentheses


            parentheses = parentheses_function(s)

            data = []
            for i, (open_idx, close_idx) in enumerate(parentheses):
                all_parentheses = s[open_idx:close_idx + 1]
                first_pair_part = all_parentheses

                node = s[close_idx + 1:close_idx + 5]
                last_idx = close_idx + 4
                while last_idx + 1 < len(s) and s[last_idx + 1].isdigit():
                    node += s[last_idx + 1]
                    last_idx += 1

                num_letters = ""
                j = last_idx
                while j < len(s) and s[j].isdigit():  # or s[j] == ".":
                    num_letters += s[j]
                    j += 1
                s2 = s[j:]

                if s2.find(")") < s2.find("(") or s2.find("(") == -1:
                    end_node = s2[s2.find(")") + 1:s2.find(")") + 5]
                    for char in s2[s2.find(")") + 5:]:
                        if char.isdigit():
                            end_node += char
                        else:
                            break
                    data.append([first_pair_part, node, end_node])

                elif s2.find("(") < s2.find(")") and s2.find("(") >= 0:
                    stack = []
                    parentheses = []

                    for i, char in enumerate(s2):
                        if char == '(':
                            stack.append(i)
                        elif char == ')':
                            if stack:
                                open_idx = stack.pop()
                                parentheses.append((open_idx, i))

                    # Find the smallest open_idx
                    smallest_open_idx = min([p[0] for p in parentheses])

                    # Print the smallest pair
                    for open_idx, close_idx in parentheses:
                        if open_idx == smallest_open_idx:
                            closed_parenthesis_pair = s2[open_idx:close_idx + 1]

                            # Print remaining chars in s
                            remaining_chars = s2[close_idx + 1:]
                            # print(remaining_chars)
                            idx = remaining_chars.find(')') + 1
                            output = remaining_chars[idx:idx + 4]
                            for char in remaining_chars[idx + 4:]:
                                if char.isdigit():
                                    output += char
                                else:
                                    break
                            end_node = output
                            data.append([first_pair_part, node, end_node])
                            # print(end_node)
            pd.set_option('display.max_colwidth', None)
            df1 = pd.DataFrame(data, columns=['leafs_one_side', 'NodeA', 'NodeB'])[:-1]

            df1['NodeA-NodeB'] = df1['NodeA'] + '-' + df1['NodeB']

            b = treefile_content

            boot_data = []
            parentheses = parentheses_function(b)
            for i, (open_idx, close_idx) in enumerate(parentheses):
                parentheses_pair = b[open_idx:close_idx + 1]
                last_idx = close_idx
                boot = ""
                while last_idx + 1 < len(b) and b[last_idx + 1].isdigit():
                    boot += b[last_idx + 1]
                    last_idx += 1
                boot_data.append([boot])
            df_boot = pd.DataFrame(boot_data, columns=['boots_value'])[:-1]

            df1 = df1.merge(df_boot, left_index=True, right_index=True)

            df1_indices = {row['NodeA-NodeB']: i for i, row in df1.iterrows()}
            df2['df1_index'] = df2['Node_A-Node_B'].map(df1_indices)
            df2_sorted = df2.sort_values(by='df1_index')
            merged_df = df2_sorted.merge(df1, left_on='df1_index', right_index=True, how='left')

            merged_df = merged_df[['delta', 'c_s', 'p_value', 'result_test', 'boots_value', 'Node_A', 'Node_B']]
            merged_df.fillna("", inplace=True)


        else:
            merged_df = df2[['delta', 'c_s', 'p_value', 'result_test', 'Node_A', 'Node_B']]

        newick = s

        ts = TreeStyle()
        ts.show_leaf_name = True
        t = Tree(newick, format=1)

        if args.midpoint:
            t.set_outgroup(t.get_midpoint_outgroup())


        def custom_layout(node):
            if node.is_leaf():
                node.dist = 0.1
            else:
                node.dist = 0.1

        if args.compact:
            ts.layout_fn = custom_layout

        canvas_width = 1200
        canvas_height = 900

        ts.scale = canvas_width / t.get_distance(t.get_leaves()[0], t.get_leaves()[-1])

        if args.cycle:
            ts.mode = "c"

        labeled_branches = set()


        if matches:
            def add_boot_label(node, boot, fgcolor):
                node.add_face(TextFace(f"{boot:}", fsize=35, fgcolor=fgcolor), column=0, position="branch-top")


        def add_delta_label(node, delta, fgcolor):
            node.add_face(TextFace(f"{delta:.4f}", fsize=35, fgcolor=fgcolor), column=0, position="branch-top")


        def add_p_label(node, p_value, fgcolor):
            node.add_face(TextFace(f"{p_value:.2e}", fsize=35, fgcolor=fgcolor), column=0, position="branch-top")


        def add_saturated_label(node, satur, fgcolor):
            node.add_face(TextFace(f"{branch_status}", fsize=35, fgcolor=fgcolor), column=0, position="branch-top")


        for _, row in merged_df.iterrows():
            node_b = row['Node_A']
            delta = float(row['delta'])
            branch_status = row['result_test']

            if matches:
                boots_value = row['boots_value']
                boot = '' if boots_value == '' else int(boots_value)  # Replace empty strings with 0


            p_value = float(row['p_value'])

            node_b_obj = t.search_nodes(name=node_b)
            if node_b_obj:
                node_b_obj = node_b_obj[0]
                if node_b_obj.name not in labeled_branches:

                    if matches:
                        add_boot_label(node_b_obj, boot, "black")


                    add_delta_label(node_b_obj, delta, "green")
                    add_p_label(node_b_obj, p_value, fgcolor="blue")

                    labeled_branches.add(node_b_obj.name)

                if branch_status == 'Saturated':
                    node_b_obj.dist += 0  # Increase the branch length
                    node_style = NodeStyle()
                    node_style["hz_line_color"] = "red"
                    node_style["hz_line_width"] = 4
                    node_style["hz_line_type"] = 2
                    node_b_obj.set_style(node_style)

                    add_saturated_label(node_b_obj, branch_status, fgcolor="red")

        for leaf in t.iter_leaves():
            delta_label = merged_df[merged_df['Node_B'] == leaf.name]['delta'].values

        for node in t.traverse():
            node_style = NodeStyle()
            node_style["size"] = 10
            node_style["fgcolor"] = "black"
            node_style["hz_line_width"] = 4
            node_style["vt_line_width"] = 4

            branch_status = merged_df[merged_df['Node_A'] == node.name]['result_test'].values
            if len(branch_status) > 0 and branch_status[0] == 'Saturated':
                node_style["hz_line_color"] = "red"
                node_style["hz_line_width"] = 4
                node_style["hz_line_type"] = 2
                node_style["size"] = 10
                node_style["fgcolor"] = "red"

            node.set_style(node_style)

            if node.is_leaf():
                node.add_face(TextFace(node.name, fsize=40), column=0, position="branch-right")
                node.name = ""

            if not node.is_leaf():
                node.add_face(TextFace(node.name, fsize=40), column=0, position="branch-right")

            branch_length = node.dist
            node.add_face(TextFace(f" {branch_length:.4f}", fsize=35, fgcolor="red"), column=1,
                          position="branch-bottom")

        # t.show(tree_style=ts)
        base_file_name = os.path.splitext(csv_file_name)[0]

        suffixes = []
        if args.midpoint:
            suffixes.append("midpoint")
        if args.cycle:
            suffixes.append("cycle")
        if args.compact:
            suffixes.append("compact")

        suffix = "_".join(suffixes)

        pdf_suffix = ".pdf"
        svg_suffix = ".svg"

        if suffix:
            pdf_suffix = f"_{suffix}{pdf_suffix}"
            svg_suffix = f"_{suffix}{svg_suffix}"

        pdf_file_name = base_file_name + pdf_suffix
        pdf_file_path = os.path.join(image_directory, pdf_file_name)

        dpi = 600  # Adjust this value as needed

        t.render(pdf_file_path, tree_style=ts, dpi=dpi)
        svg_file_name = base_file_name + svg_suffix
        svg_file_path = os.path.join(image_directory, svg_file_name)
        t.render(svg_file_path, tree_style=ts, dpi=dpi)


else:
    print("No matching .treefile found.")


















