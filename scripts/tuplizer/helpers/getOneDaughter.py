def getOneDaughter(mesons, GenPart_genPartIdxMother_, nGenPart_):
    """
    Finds one daughter particle for each meson in the provided list.

    Parameters:
        mesons (array-like): List of meson indices.
        GenPart_genPartIdxMother_ (array-like): Array indicating the mother index of each generated particle.
        nGenPart_ (int): Total number of generated particles.

    Returns:
        list: A list containing the index of one daughter particle for each meson.
              If no daughter is found, the value -1 is assigned.
    """
    oneDaughter = []      # one index of a daughter of the meson aligned
    for mes in mesons:
        foundDaughter = -1
        for gp in range(nGenPart_):
            if (GenPart_genPartIdxMother_[gp] == mes):
                oneDaughter.append(gp)
                foundDaughter=gp
                # fill one daughter per meson
                break
        if foundDaughter==-1:
            oneDaughter.append(-1)
    assert len(mesons)==len(oneDaughter)
    return oneDaughter