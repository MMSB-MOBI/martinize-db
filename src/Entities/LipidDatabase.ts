import { Lipid } from "./entities";
import AbstractDatabase from "./AbstractDatabase";
import RadiusDatabase from "./RadiusDatabase";

export default class LipidDatabase extends AbstractDatabase<Lipid> {
  async getAndThrowIfMissing(keys: string[], force_field: string) {
    const prefix = this.getPrefix(force_field);

    const unique = new Set(keys.map(e => prefix + e));
    const items = await this.bulkGet([...unique]);

    const fetched = new Set(items.map(i => i.id));

    for (const item of unique) {
      if (!fetched.has(item)) {
        throw new Error("Lipid " + item + " is missing.");
      }
    }

    return items;
  }

  protected getPrefix(force_field: string) {
    const prefix = RadiusDatabase.FORCE_FIELD_TO_MARTINI_VERSION[force_field];
    if (!prefix) {
      throw new Error('Unsupported force field.');
    }

    return prefix + '~';
  }

  async getAvailableNames(force_field?: string) {
    const keys = await this.keys();

    if (force_field) {
      const prefix = this.getPrefix(force_field);

      const names = keys.filter(e => e.startsWith(prefix)).map(e => e.split('~')[1]);

      return [...new Set(names)];
    }

    const names = keys.map(e => e.split('~')[1]);

    return [...new Set(names)];
  }

  async add(lipid: Lipid, force_field: string) {
    lipid.id = this.getPrefix(force_field) + lipid.name;

    const exists = await this.exists(lipid.id);

    if (!exists)
      await this.save(lipid);
  }

  async removeByName(name: string, force_field: string) {
    const id = this.getPrefix(force_field) + name;
    return this.delete(await this.get(id));
  }
}
