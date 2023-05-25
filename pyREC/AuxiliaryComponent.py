class AuxiliaryComponent:
    def __init__(self, id, tech, replacement_year, replacement_cost):
        """
        :param id: identification code
        :param tech:'storage', 'pump','compressor
        :param replacement_year
        :param replacement_cost: initial cost [€]
        """
        self.id = id
        self.tech = tech
        self.replacement_year = replacement_year
        self.replacement_cost = replacement_cost
