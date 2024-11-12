def parse_ac_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = {
        'CHARGE': None,
        'Formula': None,
        'Atoms': [],
        'Bonds': []
    }

    for line in lines:
        parts = line.split()
        if not parts:
            continue

        if parts[0] == 'ATOM':
            atom_info = {
                'Type': parts[2],
                'Residue': parts[9],
                'Charge': float(parts[8]),
                'Index': int(parts[1])
            }
            data['Atoms'].append(atom_info)
        elif parts[0] == 'BOND':
            bond_info = {
                'Atom1': parts[5],
                'Atom2': parts[6]
            }
            data['Bonds'].append(bond_info)

    return data

def format_rtp(data):
    formatted = []

    # Add header
    formatted.append("[ C16 ]")
    formatted.append(" [ atoms ]")
    for atom in data['Atoms']:
        formatted.append(f"  {atom['Type']:>4} {atom['Residue']:>4}       {atom['Charge']:>12.5f}  {atom['Index']:>4}")

    # Add bonds
    formatted.append("\n [ bonds ]")
    for bond in data['Bonds']:
        formatted.append(f"  {bond['Atom1']:>4}  {bond['Atom2']:>4}")

    # Add ending and impropers
    formatted.append("    -C     N")
    formatted.append("\n [ impropers ]")
    formatted.append("    -C    CA     N     H")
    formatted.append("    CA    +N     C     O")

    return '\n'.join(formatted)

def save_rtp_file(data, output_path):
    formatted_data = format_rtp(data)
    with open(output_path, 'w') as file:
        file.write(formatted_data)

if __name__ == "__main__":
    input_file = '/Users/annettethompson/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/MD/KS_Mutants/FabF_WT_dimer_C16/Amber C16/C16.ac'
    output_file = '/Users/annettethompson/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/MD/KS_Mutants/FabF_WT_dimer_C16/Amber C16/C16_rtp.txt'
    data = parse_ac_file(input_file)
    save_rtp_file(data, output_file)
    print(f"Data saved to {output_file}")
